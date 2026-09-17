import time
import uuid
import json
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Any, List, Optional, Generator

from backend.services.batch import BatchScreeningService

class BatchScreeningManager:
    """
    Manages asynchronous virtual screening jobs in the background using ThreadPoolExecutor.
    Supports real-time status polling, Server-Sent Events (SSE) streaming, and cancellation.
    """
    _instance = None
    _lock = threading.Lock()

    def __init__(self, max_workers: int = 2):
        self.executor = ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="batch_docking_worker")
        self.jobs: Dict[str, Dict[str, Any]] = {}

    @classmethod
    def get_instance(cls) -> "BatchScreeningManager":
        with cls._lock:
            if cls._instance is None:
                cls._instance = BatchScreeningManager()
            return cls._instance

    def start_batch_job(
        self,
        receptor_pdbqt: str,
        receptor_pdb: str,
        pocket_center: Dict[str, float],
        pocket_size: Dict[str, float],
        ligand_list: List[Dict[str, str]],
        exhaustiveness: int = 4
    ) -> str:
        job_id = f"batch_{int(time.time())}_{uuid.uuid4().hex[:8]}"
        job_data = {
            "job_id": job_id,
            "status": "queued",  # queued, running, completed, failed, cancelled
            "total": len(ligand_list),
            "completed": 0,
            "percent": 0.0,
            "current_compound": "",
            "leaderboard": [],
            "error": None,
            "cancel_requested": False,
            "created_at": time.time(),
            "updated_at": time.time()
        }

        with self._lock:
            self.jobs[job_id] = job_data

        self.executor.submit(
            self._run_job_worker,
            job_id,
            receptor_pdbqt,
            receptor_pdb,
            pocket_center,
            pocket_size,
            ligand_list,
            exhaustiveness
        )

        return job_id

    def get_job(self, job_id: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            job = self.jobs.get(job_id)
            if job:
                # Return shallow copy
                return dict(job)
            return None

    def cancel_job(self, job_id: str) -> bool:
        with self._lock:
            job = self.jobs.get(job_id)
            if job and job["status"] in ("queued", "running"):
                job["cancel_requested"] = True
                job["status"] = "cancelled"
                job["updated_at"] = time.time()
                return True
            return False

    def stream_job_events(self, job_id: str) -> Generator[str, None, None]:
        """Yield Server-Sent Events (SSE) formatted data for real-time browser updates."""
        last_completed = -1
        while True:
            job = self.get_job(job_id)
            if not job:
                yield f"data: {json.dumps({'error': 'Job not found', 'status': 'failed'})}\n\n"
                break

            status = job.get("status")
            completed = job.get("completed", 0)

            # Send update on state change or progress increment
            if completed != last_completed or status in ("completed", "failed", "cancelled"):
                last_completed = completed
                event_data = {
                    "job_id": job_id,
                    "status": status,
                    "total": job["total"],
                    "completed": completed,
                    "percent": job["percent"],
                    "current_compound": job["current_compound"],
                    "leaderboard": job["leaderboard"],
                    "error": job["error"]
                }
                yield f"data: {json.dumps(event_data)}\n\n"

            if status in ("completed", "failed", "cancelled"):
                break

            time.sleep(0.4)

    def _run_job_worker(
        self,
        job_id: str,
        receptor_pdbqt: str,
        receptor_pdb: str,
        pocket_center: Dict[str, float],
        pocket_size: Dict[str, float],
        ligand_list: List[Dict[str, str]],
        exhaustiveness: int
    ) -> None:
        with self._lock:
            job = self.jobs.get(job_id)
            if not job or job["cancel_requested"]:
                return
            job["status"] = "running"
            job["updated_at"] = time.time()

        def progress_callback(completed_count: int, total_count: int, current_name: str, partial_leaderboard: List[Dict[str, Any]]):
            with self._lock:
                j = self.jobs.get(job_id)
                if j:
                    j["completed"] = completed_count
                    j["total"] = total_count
                    j["percent"] = round((completed_count / max(1, total_count)) * 100.0, 1)
                    j["current_compound"] = current_name
                    j["leaderboard"] = partial_leaderboard
                    j["updated_at"] = time.time()

        def is_cancelled() -> bool:
            with self._lock:
                j = self.jobs.get(job_id)
                return j.get("cancel_requested", False) if j else True

        try:
            leaderboard = BatchScreeningService.run_batch(
                receptor_pdbqt=receptor_pdbqt,
                receptor_pdb=receptor_pdb,
                pocket_center=pocket_center,
                pocket_size=pocket_size,
                ligand_list=ligand_list,
                exhaustiveness=exhaustiveness,
                progress_callback=progress_callback,
                is_cancelled_fn=is_cancelled
            )
            with self._lock:
                j = self.jobs.get(job_id)
                if j and not j["cancel_requested"]:
                    j["status"] = "completed"
                    j["completed"] = len(ligand_list)
                    j["percent"] = 100.0
                    j["leaderboard"] = leaderboard
                    j["updated_at"] = time.time()
        except Exception as e:
            with self._lock:
                j = self.jobs.get(job_id)
                if j:
                    j["status"] = "failed"
                    j["error"] = str(e)
                    j["updated_at"] = time.time()
