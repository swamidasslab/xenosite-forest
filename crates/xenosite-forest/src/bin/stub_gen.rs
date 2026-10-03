//! Write (default) or verify (`--check`) `src/xenosite/forest/_rust.pyi`.
//!
//! `make stubs` / `make check-stubs` (`stubs` feature). The unit test
//! `checked_in_stub_is_current` is the same drift gate.

use std::process::ExitCode;

fn main() -> ExitCode {
    let check = std::env::args().any(|a| a == "--check");
    let want = match xenosite_forest::python_stubs::render_stub() {
        Ok(text) => text,
        Err(e) => {
            eprintln!("stub_gen: {e:#}");
            return ExitCode::FAILURE;
        }
    };
    let path = xenosite_forest::python_stubs::stub_path();
    let rel = xenosite_forest::python_stubs::STUB_PATH;
    if check {
        let have = std::fs::read_to_string(&path).unwrap_or_default();
        if have != want {
            eprintln!("stub_gen: {rel} is stale — run `make stubs` and commit the result");
            return ExitCode::FAILURE;
        }
        println!("stub_gen: {rel} is current");
        return ExitCode::SUCCESS;
    }
    if let Err(e) = std::fs::write(&path, want) {
        eprintln!("stub_gen: write {rel}: {e}");
        return ExitCode::FAILURE;
    }
    println!("stub_gen: wrote {rel}");
    ExitCode::SUCCESS
}
