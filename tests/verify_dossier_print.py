import asyncio
from pathlib import Path
from playwright.async_api import async_playwright

async def run():
    async with async_playwright() as p:
        browser = await p.chromium.launch()
        page = await browser.new_page()
        await page.goto("http://127.0.0.1:5000", wait_until="networkidle")
        await page.wait_for_selector("#btn-export-dossier", state="attached")
        
        # Populate state with sample docking session (3ERT / Tamoxifen like user's screenshot)
        js_code = """
        () => {
            const app = window.bindoraApp;
            app.state.receptor = {
                pdb_id: '3ERT',
                title: 'HUMAN ESTROGEN RECEPTOR ALPHA LIGAND-BINDING DOMAIN IN COMPLEX WITH 4-HYDROXYTAMOXIFEN',
                atom_count: 1950,
                chains: ['A'],
                resolution: 1.90
            };
            app.state.ligand = {
                name: 'Tamoxifen',
                formula: 'C26H29NO',
                canonical_smiles: 'CCC(=C(C1=CC=CC=C1)C2=CC=C(C=C2)OCCN(C)C)C3=CC=CC=C3',
                heavy_atom_count: 28,
                weight: 371.51,
                adme: {
                    physicochemical: {
                        molecular_weight: { value: 371.51 },
                        logp: { value: 6.30 },
                        h_bond_donors: { value: 0 },
                        h_bond_acceptors: { value: 2 },
                        rotatable_bonds: { value: 8 },
                        tpsa: { value: 12.47 }
                    },
                    drug_likeness: {
                        lipinski: { status: 'Borderline', violations_count: 1, violations: ['LogP > 5.0'] }
                    },
                    pharmacokinetics: {
                        gi_absorption: { level: 'High' },
                        bbb_permeant: { is_permeant: true },
                        cyp450_inhibition: { alerts: [{ cyp: 'CYP2D6' }, { cyp: 'CYP3A4' }] }
                    }
                }
            };
            app.state.docking = {
                top_pose: {
                    mode: 1,
                    affinity_kcal: -9.75,
                    vinardo_affinity_kcal: -7.21,
                    rmsd_lb: 0.0,
                    rmsd_ub: 0.0
                },
                poses: [
                    { mode: 1, affinity_kcal: -9.75, vinardo_affinity_kcal: -7.21, rmsd_lb: 0.0, rmsd_ub: 0.0 },
                    { mode: 2, affinity_kcal: -9.42, vinardo_affinity_kcal: -6.98, rmsd_lb: 1.24, rmsd_ub: 1.85 },
                    { mode: 3, affinity_kcal: -9.10, vinardo_affinity_kcal: -6.74, rmsd_lb: 1.89, rmsd_ub: 2.31 },
                    { mode: 4, affinity_kcal: -8.85, vinardo_affinity_kcal: -6.50, rmsd_lb: 2.15, rmsd_ub: 2.94 },
                    { mode: 5, affinity_kcal: -8.60, vinardo_affinity_kcal: -6.32, rmsd_lb: 2.48, rmsd_ub: 3.42 },
                    { mode: 6, affinity_kcal: -8.31, vinardo_affinity_kcal: -6.11, rmsd_lb: 2.91, rmsd_ub: 3.90 },
                    { mode: 7, affinity_kcal: -8.05, vinardo_affinity_kcal: -5.92, rmsd_lb: 3.12, rmsd_ub: 4.15 },
                    { mode: 8, affinity_kcal: -7.82, vinardo_affinity_kcal: -5.74, rmsd_lb: 3.45, rmsd_ub: 4.52 },
                    { mode: 9, affinity_kcal: -7.50, vinardo_affinity_kcal: -5.50, rmsd_lb: 3.80, rmsd_ub: 4.95 }
                ],
                thermodynamics: {
                    theoretical_kd_nm: 71.3,
                    ligand_efficiency: { value: 0.348 },
                    size_independent_le: { value: 1.84 }
                },
                interactions: {
                    total_hbond_count: 1,
                    total_hydrophobic_count: 4,
                    hydrogen_bonds: [
                        { residue: 'GLU 353:A', res_name: 'GLU', res_num: '353', chain: 'A', receptor_atom: 'OE1', ligand_atom: 'O1', distance: 2.74, type: 'Hydrogen Bond' }
                    ],
                    hydrophobic_contacts: [
                        { residue: 'LEU 387:A', res_name: 'LEU', res_num: '387', chain: 'A', distance: 3.45 },
                        { residue: 'MET 388:A', res_name: 'MET', res_num: '388', chain: 'A', distance: 3.62 },
                        { residue: 'ILE 424:A', res_name: 'ILE', res_num: '424', chain: 'A', distance: 3.78 },
                        { residue: 'PHE 404:A', res_name: 'PHE', res_num: '404', chain: 'A', distance: 3.82 }
                    ]
                },
                grid: {
                    center: { x: 30.1, y: -1.8, z: 24.3 },
                    size: { x: 22.0, y: 22.0, z: 22.0 }
                },
                exhaustiveness: 16
            };
            app.state.redockingValidation = {
                pdb_id: '3ERT',
                rmsd_angstroms: 0.91,
                is_validated: true,
                benchmark_status: 'Pass (Publication Grade)'
            };
            app.updateDossierView();
            // Switch to tab-dossier
            document.querySelectorAll('.tab-pane').forEach(el => el.classList.add('hidden'));
            document.getElementById('tab-dossier').classList.remove('hidden');
        }
        """
        await page.evaluate(js_code)
        
        # Emulate print media
        await page.emulate_media(media="print")
        
        # Take print screenshot
        screenshot_path = "tests/dossier_print_preview.png"
        await page.screenshot(path=screenshot_path, full_page=True)
        print("Screenshot saved:", screenshot_path)
        
        # Also generate PDF
        pdf_path = "tests/dossier_export.pdf"
        await page.pdf(path=pdf_path, format="A4", print_background=True)
        print("PDF saved:", pdf_path)
        
        await browser.close()

if __name__ == "__main__":
    asyncio.run(run())

