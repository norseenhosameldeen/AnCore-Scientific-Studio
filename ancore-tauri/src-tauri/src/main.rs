#![cfg_attr(
    all(not(debug_assertions), target_os = "macos"),
    windows_subsystem = "windows"
)]

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Serialize, Deserialize)]
struct TargetBundle {
    pdb_id: String,
    target_name: String,
    ligand_resname: String,
    clinical_drug: String,
    is_custom: bool,
    has_all_stages: bool,
    stage1: Value,
    stage2: Value,
    stage3: Value,
    stage4: Value,
    stage5: Value,
    stage6: Value,
    receptor_pdb: String,
    crystal_ligand_pdb: String,
}

fn read_json_safe(p: &Path) -> Value {
    if p.exists() {
        if let Ok(txt) = fs::read_to_string(p) {
            if let Ok(v) = serde_json::from_str::<Value>(&txt) {
                return v;
            }
        }
    }
    json!({})
}

fn read_text_safe(p_str: Option<&str>) -> String {
    if let Some(s) = p_str {
        let p = Path::new(s);
        if p.exists() {
            if let Ok(txt) = fs::read_to_string(p) {
                return txt;
            }
        }
    }
    String::new()
}

#[tauri::command]
fn load_workspace(project_dir: String) -> Result<Value, String> {
    let proj = PathBuf::from(&project_dir);
    let bench_root = proj.join("benchmarks");

    let default_targets = vec![
        ("1IEP", "ABL1 Tyrosine Kinase", "STI", "Imatinib"),
        ("3ERT", "Estrogen Receptor alpha", "OHT", "4-Hydroxytamoxifen"),
        ("1UYG", "HSP90-alpha Chaperone", "PU3", "Purine PU3"),
        ("1M17", "EGFR Tyrosine Kinase", "AQ4", "Erlotinib"),
        ("1STP", "Streptavidin", "BTN", "D-Biotin"),
    ];

    let benchmark_ids: Vec<String> = default_targets.iter().map(|t| t.0.to_string()).collect();
    let mut known_ids: Vec<String> = benchmark_ids.clone();

    if bench_root.exists() {
        if let Ok(entries) = fs::read_dir(&bench_root) {
            for entry in entries.flatten() {
                if let Ok(ft) = entry.file_type() {
                    if ft.is_dir() {
                        let name = entry.file_name().to_string_lossy().to_uppercase();
                        if name.len() >= 3 && name.len() <= 8 && !name.starts_with('_') && !known_ids.contains(&name) {
                            known_ids.push(name);
                        }
                    }
                }
            }
        }
    }

    let mut targets_map = serde_json::Map::new();
    let mut target_order: Vec<String> = Vec::new();
    let mut custom_order: Vec<String> = Vec::new();

    for pdb_id in known_ids {
        let tdir = bench_root.join(&pdb_id);
        if !tdir.exists() {
            continue;
        }

        let s1 = read_json_safe(&tdir.join("stage1_manifest.json"));
        let s2 = read_json_safe(&tdir.join("stage2_manifest.json"));
        let s3 = read_json_safe(&tdir.join("stage3_manifest.json"));
        let s4 = read_json_safe(&tdir.join("stage4_manifest.json"));
        let s5 = read_json_safe(&tdir.join("stage5_manifest.json"));
        let s6 = read_json_safe(&tdir.join("stage6_manifest.json"));

        let has_all_stages = tdir.join("stage1_manifest.json").exists()
            && tdir.join("stage6_manifest.json").exists();

        let def_meta = default_targets.iter().find(|x| x.0 == pdb_id);
        let is_custom = def_meta.is_none();

        let target_name = s1
            .get("target_name")
            .and_then(|v| v.as_str())
            .or_else(|| def_meta.map(|x| x.1))
            .unwrap_or(&format!("Target {}", pdb_id))
            .to_string();

        let ligand_resname = s1
            .get("ligand_id")
            .and_then(|v| v.as_str())
            .or_else(|| def_meta.map(|x| x.2))
            .unwrap_or("LIG")
            .to_string();

        let clinical_drug = s1
            .get("ligand_name")
            .and_then(|v| v.as_str())
            .or_else(|| def_meta.map(|x| x.3))
            .unwrap_or(&ligand_resname)
            .to_string();

        let rec_path = s1
            .get("files")
            .and_then(|f| f.get("receptor_clean_pdb"))
            .and_then(|v| v.as_str());
        let cryst_path = s1
            .get("files")
            .and_then(|f| f.get("crystal_ligand_pdb"))
            .and_then(|v| v.as_str());

        let mut receptor_pdb = read_text_safe(rec_path);
        if receptor_pdb.is_empty() {
            receptor_pdb = read_text_safe(tdir.join("receptor_clean.pdb").to_str());
        }
        let mut crystal_ligand_pdb = read_text_safe(cryst_path);
        if crystal_ligand_pdb.is_empty() {
            crystal_ligand_pdb = read_text_safe(tdir.join("crystal_ligand.pdb").to_str());
        }

        let bundle = TargetBundle {
            pdb_id: pdb_id.clone(),
            target_name,
            ligand_resname,
            clinical_drug,
            is_custom,
            has_all_stages,
            stage1: s1,
            stage2: s2,
            stage3: s3,
            stage4: s4,
            stage5: s5,
            stage6: s6,
            receptor_pdb,
            crystal_ligand_pdb,
        };

        if let Ok(val) = serde_json::to_value(bundle) {
            targets_map.insert(pdb_id.clone(), val);
            target_order.push(pdb_id.clone());
            if is_custom {
                custom_order.push(pdb_id);
            }
        }
    }

    Ok(json!({
        "project_dir": project_dir,
        "target_order": target_order,
        "benchmark_order": benchmark_ids,
        "custom_order": custom_order,
        "targets": targets_map
    }))
}

#[tauri::command]
fn run_python_pipeline(
    project_dir: String,
    python_bin: String,
    pdb_id: String,
    target_name: String,
    ligand_code: String,
    custom_smiles: String,
    local_pdb_content: Option<String>,
    sdf_content: Option<String>,
    stage_only: Option<i32>,
    reset_downstream: Option<bool>,
    ph: Option<f64>,
    box_padding: Option<f64>,
    exhaustiveness: Option<i32>,
    temp_k: Option<f64>,
    rmsd_gate: Option<f64>,
) -> Result<String, String> {
    let proj = PathBuf::from(&project_dir);
    let runner_script = if proj.join("src").join("ancore_custom_runner.py").exists() {
        proj.join("src").join("ancore_custom_runner.py")
    } else {
        proj.join("ancore_custom_runner.py")
    };

    let code = pdb_id.trim().to_uppercase();
    let tdir = proj.join("benchmarks").join(&code);
    let _ = fs::create_dir_all(&tdir);

    if let Some(pdb_txt) = local_pdb_content {
        if !pdb_txt.trim().is_empty() {
            let raw_file = tdir.join(format!("{}_raw.pdb", code));
            fs::write(&raw_file, pdb_txt)
                .map_err(|e| format!("Failed to save local PDB file: {}", e))?;
        }
    }

    if let Some(sdf_txt) = sdf_content {
        if !sdf_txt.trim().is_empty() {
            let sdf_file = tdir.join("uploaded_ligands.sdf");
            fs::write(&sdf_file, sdf_txt)
                .map_err(|e| format!("Failed to save uploaded SDF file: {}", e))?;
        }
    }

    let py_cmd = if Path::new(&python_bin).exists() {
        python_bin
    } else {
        "python3".to_string()
    };

    let mut cmd = Command::new(&py_cmd);
    cmd.arg(&runner_script)
        .arg("--project-dir")
        .arg(&project_dir)
        .arg("--pdb-id")
        .arg(&pdb_id)
        .arg("--target-name")
        .arg(&target_name)
        .arg("--ligand-code")
        .arg(&ligand_code)
        .arg("--custom-smiles")
        .arg(&custom_smiles)
        .arg("--stage")
        .arg(stage_only.unwrap_or(1).to_string())
        .arg("--reset-downstream")
        .arg(if reset_downstream.unwrap_or(false) { "1" } else { "0" })
        .arg("--ph")
        .arg(ph.unwrap_or(7.4).to_string())
        .arg("--padding")
        .arg(box_padding.unwrap_or(4.0).to_string())
        .arg("--exhaustiveness")
        .arg(exhaustiveness.unwrap_or(8).to_string())
        .arg("--temp-k")
        .arg(temp_k.unwrap_or(300.0).to_string())
        .arg("--rmsd-gate")
        .arg(rmsd_gate.unwrap_or(2.0).to_string())
        .current_dir(&proj);

    let output = cmd
        .output()
        .map_err(|e| format!("Failed to launch scientific engine: {}", e))?;

    if output.status.success() {
        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    } else {
        let stderr_txt = String::from_utf8_lossy(&output.stderr).to_string();
        let stdout_txt = String::from_utf8_lossy(&output.stdout).to_string();
        Err(format!("{}\n{}", stderr_txt, stdout_txt))
    }
}

#[tauri::command]
fn delete_custom_target(project_dir: String, pdb_id: String) -> Result<String, String> {
    let proj = PathBuf::from(&project_dir);
    let code = pdb_id.trim().to_uppercase();
    let default_set = ["1IEP", "3ERT", "1UYG", "1M17", "1STP"];
    if default_set.contains(&code.as_str()) {
        return Err("Cannot delete a built-in benchmark reference target.".to_string());
    }
    let tdir = proj.join("benchmarks").join(&code);
    if tdir.exists() {
        fs::remove_dir_all(&tdir).map_err(|e| format!("Failed to remove {}: {}", code, e))?;
    }
    Ok(format!("Deleted {}", code))
}

fn main() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![
            load_workspace,
            run_python_pipeline,
            delete_custom_target
        ])
        .run(tauri::generate_context!())
        .expect("error while running AnCore Studio Tauri application");
}
