// Chart.js visualizations for ADME, Lipinski Radar, and Physicochemical properties

class ADMECharts {
  static radarChartInstance = null;

  static renderLipinskiRadar(canvasId, admeData) {
    const canvas = document.getElementById(canvasId);
    if (!canvas || !window.Chart) return;

    if (this.radarChartInstance) {
      this.radarChartInstance.destroy();
      this.radarChartInstance = null;
    }

    if (!admeData || !admeData.physicochemical) return;

    const p = admeData.physicochemical;
    const mw = p.molecular_weight?.value || 0;
    const logp = p.logp?.value || 0;
    const hbd = p.hbd?.value || 0;
    const hba = p.hba?.value || 0;
    const rotb = p.rotatable_bonds?.value || 0;
    const tpsa = p.tpsa?.value || 0;
    const fsp3 = p.fsp3?.value || 0;

    // Normalize each descriptor onto a 0-100 scale where 100 represents max acceptable boundary
    // MW: 0-500 Da
    const normMW = Math.min(100, Math.max(0, (mw / 500) * 100));
    // LogP: -2 to 5 -> 0-100
    const normLogP = Math.min(100, Math.max(0, ((logp + 2) / 7) * 100));
    // HBD: 0-5 -> 0-100
    const normHBD = Math.min(100, Math.max(0, (hbd / 5) * 100));
    // HBA: 0-10 -> 0-100
    const normHBA = Math.min(100, Math.max(0, (hba / 10) * 100));
    // RotB: 0-10 -> 0-100
    const normRotB = Math.min(100, Math.max(0, (rotb / 10) * 100));
    // TPSA: 0-140 -> 0-100
    const normTPSA = Math.min(100, Math.max(0, (tpsa / 140) * 100));

    const labels = [
      `MW (${mw} Da)`,
      `LogP (${logp})`,
      `HBD (${hbd})`,
      `HBA (${hba})`,
      `RotB (${rotb})`,
      `TPSA (${tpsa} Å²)`
    ];

    const ctx = canvas.getContext("2d");
    this.radarChartInstance = new Chart(ctx, {
      type: "radar",
      data: {
        labels: labels,
        datasets: [
          {
            label: "Drug Candidate",
            data: [normMW, normLogP, normHBD, normHBA, normRotB, normTPSA],
            backgroundColor: "rgba(56, 189, 248, 0.35)", // Cyan
            borderColor: "#38bdf8",
            borderWidth: 2,
            pointBackgroundColor: "#0284c7",
            pointBorderColor: "#ffffff",
            pointRadius: 4
          },
          {
            label: "Lipinski Rule Boundary (Threshold)",
            data: [100, 100, 100, 100, 100, 100],
            backgroundColor: "rgba(16, 185, 129, 0.08)", // Emerald
            borderColor: "#10b981",
            borderDash: [4, 4],
            borderWidth: 1.5,
            pointRadius: 0
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          r: {
            angleLines: { color: "rgba(255, 255, 255, 0.1)" },
            grid: { color: "rgba(255, 255, 255, 0.08)" },
            pointLabels: {
              color: "#cbd5e1",
              font: { size: 11, family: "sans-serif" }
            },
            ticks: {
              display: false,
              min: 0,
              max: 120
            }
          }
        },
        plugins: {
          legend: {
            labels: {
              color: "#94a3b8",
              font: { size: 11 }
            }
          },
          tooltip: {
            backgroundColor: "rgba(15, 23, 42, 0.9)",
            titleColor: "#38bdf8",
            bodyColor: "#ffffff"
          }
        }
      }
    });
  }
}

window.ADMECharts = ADMECharts;
