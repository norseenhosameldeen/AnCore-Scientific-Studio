// ======================================================================================
// ANCORE SCIENTIFIC STUDIO — NATIVE RUST / TAURI DESKTOP BACKEND (`src-tauri/src/main.rs`)
// ======================================================================================

#![cfg_attr(
    all(not(debug_assertions), target_os = "windows"),
    windows_subsystem = "windows"
)]

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Serialize, Deserialize)]
pub struct Complex3DPayload {
    pub pdb_id: String,
    pub ligand_id: String,
    pub pdb_text: String,
    pub contacts_csv_text: String,
    pub pharmacophore_png_path: String,
    pub trajectory_csv_text: String,
}

fn resolve_project_root() -> PathBuf {
    if let Ok(env_root) = std::env::var("ANCORE_PROJECT_ROOT") {
        let p = PathBuf::from(env_root);
        if p.exists() {
            return p;
        }
    }
    let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    if cwd.join("workspaces").exists() || cwd.join("src").exists() {
        return cwd;
    }
    if let Some(parent) = cwd.parent() {
        if parent.join("workspaces").exists() || parent.join("src").exists() {
            return parent.to_path_buf();
        }
    }
    cwd
}

fn resolve_python_bin() -> String {
    let candidates = [
        "/Library/miniconda3/bin/python3",
        "/opt/miniconda3/bin/python3",
        "/Users/ahmedsenary/miniconda3/bin/python3",
        "/opt/homebrew/bin/python3",
        "/usr/local/bin/python3",
        "python3",
    ];
    for c in candidates.iter() {
        if Path::new(c).exists() {
            return c.to_string();
        }
    }
    "python3".to_string()
}

fn read_json_file(path: &Path) -> Value {
    match fs::read_to_string(path) {
        Ok(content) => serde_json::from_str(&content).unwrap_or_else(|_| json!(null)),
        Err(_) => json!(null),
    }
}

#[tauri::command]
fn get_suite_overview() -> Result<Value, String> {
    let root = resolve_project_root();
    let ws = root.join("workspaces");

    let s1 = read_json_file(&ws.join("stage1_suite_summary.json"));
    let s2 = read_json_file(&ws.join("stage2_suite_summary.json"));
    let s3 = read_json_file(&ws.join("stage3_suite_summary.json"));
    let s4 = read_json_file(&ws.join("stage4_suite_summary.json"));
    let s5 = read_json_file(&ws.join("stage5_suite_summary.json"));
    let s6 = read_json_file(&ws.join("stage6_suite_summary.json"));

    Ok(json!({
        "project_root": root.to_string_lossy(),
        "workspaces_dir": ws.to_string_lossy(),
        "stage1_suite": s1,
        "stage2_suite": s2,
        "stage3_suite": s3,
        "stage4_suite": s4,
        "stage5_suite": s5,
        "stage6_suite": s6,
    }))
}

#[tauri::command]
fn get_target_stage_manifests(pdb_id: String, custom_out_dir: Option<String>) -> Result<Value, String> {
    let root = resolve_project_root();
    let base_ws = match custom_out_dir {
        Some(ref d) if !d.trim().is_empty() => PathBuf::from(d.trim()),
        _ => root.join("workspaces"),
    };
    let target_dir = base_ws.join(pdb_id.trim().to_uppercase());
    if !target_dir.exists() {
        return Err(format!("Target workspace directory not found: {}", target_dir.display()));
    }

    Ok(json!({
        "pdb_id": pdb_id.to_uppercase(),
        "workspace_dir": target_dir.to_string_lossy(),
        "stage1": read_json_file(&target_dir.join("stage1_manifest.json")),
        "stage2": read_json_file(&target_dir.join("stage2_manifest.json")),
        "stage3": read_json_file(&target_dir.join("stage3_manifest.json")),
        "stage4": read_json_file(&target_dir.join("stage4_manifest.json")),
        "stage5": read_json_file(&target_dir.join("stage5_manifest.json")),
        "stage6": read_json_file(&target_dir.join("stage6_manifest.json")),
    }))
}

#[tauri::command]
fn load_3d_complex_and_contacts(pdb_id: String, ligand_id: String, custom_out_dir: Option<String>) -> Result<Complex3DPayload, String> {
    let root = resolve_project_root();
    let base_ws = match custom_out_dir {
        Some(ref d) if !d.trim().is_empty() => PathBuf::from(d.trim()),
        _ => root.join("workspaces"),
    };
    let pid = pdb_id.trim().to_uppercase();
    let lid = ligand_id.trim().to_string();
    let tdir = base_ws.join(&pid);

    let complex_pdb = if lid == "REDOCK" {
        tdir.join("02_Stage2_Redocking_Gate").join("redock_superposition_complex.pdb")
    } else {
        tdir.join("04_Stage4_Consensus_Docking")
            .join("docked_complexes_pdb")
            .join(format!("{}_complex.pdb", lid))
    };

    let pdb_text = fs::read_to_string(&complex_pdb)
        .map_err(|e| format!("Failed to read complex PDB ({}): {}", complex_pdb.display(), e))?;

    let contacts_csv = tdir
        .join("05_Stage5_Interaction_Fingerprints")
        .join("atomic_contacts_csv")
        .join(format!("{}_contacts.csv", lid));
    let contacts_csv_text = fs::read_to_string(&contacts_csv).unwrap_or_default();

    let pharm_png = tdir
        .join("05_Stage5_Interaction_Fingerprints")
        .join("pharmacophore_maps_2d")
        .join(format!("{}_pharmacophore_2d.png", lid));

    let traj_csv = tdir
        .join("06_Stage6_MD_Triage")
        .join("rmsd_trajectories_csv")
        .join(format!("{}_300K_trajectory.csv", lid));
    let trajectory_csv_text = fs::read_to_string(&traj_csv).unwrap_or_default();

    Ok(Complex3DPayload {
        pdb_id: pid,
        ligand_id: lid,
        pdb_text,
        contacts_csv_text,
        pharmacophore_png_path: pharm_png.to_string_lossy().to_string(),
        trajectory_csv_text,
    })
}

#[tauri::command]
fn run_stage_or_pipeline(script_name: String, args: Vec<String>) -> Result<Value, String> {
    let root = resolve_project_root();
    let script_path = if root.join("src").join(&script_name).exists() {
        root.join("src").join(&script_name)
    } else {
        root.join(&script_name)
    };

    if !script_path.exists() {
        return Err(format!("Python stage script not found: {}", script_path.display()));
    }

    let py_bin = resolve_python_bin();
    let output = Command::new(&py_bin)
        .arg(&script_path)
        .args(&args)
        .current_dir(&root)
        .output()
        .map_err(|e| format!("Failed to spawn {} for {}: {}", py_bin, script_name, e))?;

    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();

    Ok(json!({
        "success": output.status.success(),
        "exit_code": output.status.code(),
        "stdout": stdout,
        "stderr": stderr,
        "project_root": root.to_string_lossy(),
    }))
}

#[tauri::command]
fn export_publication_bundle(export_dir: Option<String>, pdb_id: Option<String>) -> Result<Value, String> {
    let mut args = vec!["--bundle-only".to_string()];
    if let Some(pid) = pdb_id {
        if !pid.trim().is_empty() && pid != "ALL" {
            args.push("--pdb-id".to_string());
            args.push(pid.trim().to_uppercase());
        } else {
            args.push("--run-all".to_string());
        }
    } else {
        args.push("--run-all".to_string());
    }
    if let Some(dir) = export_dir {
        if !dir.trim().is_empty() {
            args.push("--export-dir".to_string());
            args.push(dir.trim().to_string());
        }
    }
    run_stage_or_pipeline("ancore_custom_runner.py".to_string(), args)
}

fn main() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![
            get_suite_overview,
            get_target_stage_manifests,
            load_3d_complex_and_contacts,
            run_stage_or_pipeline,
            export_publication_bundle
        ])
        .run(tauri::generate_context!())
        .expect("error while running AnCore Scientific Studio Tauri application");
}