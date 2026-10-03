//! Shared local publication policy; schema writers and numerical budgets stay separate.
use std::fs;
use std::path::{Path, PathBuf};

use hdf5::File;

use super::hdf5_error;

/// The writer must release child HDF5 handles before returning. File ownership
/// stays here so every error closes the file before owned scratch is removed.
pub(crate) fn publish<T, E: From<String>>(
    output: &Path,
    inputs: &[&Path],
    overwrite: bool,
    maximum_file_bytes: Option<u64>,
    cleanup_code: &'static str,
    write: impl FnOnce(&File, &Path) -> Result<T, E>,
) -> Result<T, E> {
    protect_inputs(output, inputs).map_err(E::from)?;
    match fs::symlink_metadata(output) {
        Ok(_) if !overwrite => {
            return Err(E::from(format!(
                "output already exists: {}; pass --overwrite to replace it",
                output.display()
            )));
        }
        Err(error) if error.kind() != std::io::ErrorKind::NotFound => {
            return Err(E::from(error.to_string()));
        }
        _ => {}
    }
    let mut temporary = output.as_os_str().to_os_string();
    temporary.push(".tmp");
    let temporary = PathBuf::from(temporary);
    // Exclusive HDF5 creation is the reservation, not a check followed by truncation.
    let file = File::create_excl(&temporary).map_err(|error| E::from(hdf5_error(error)))?;
    let cleanup = TemporaryOutput {
        path: temporary,
        code: cleanup_code,
    };
    (|| {
        let file = file;
        let result = write(&file, &cleanup.path)?;
        file.flush().map_err(|error| E::from(hdf5_error(error)))?;
        file.close().map_err(|error| E::from(hdf5_error(error)))?;
        if let Some(maximum) = maximum_file_bytes {
            check_file_size(&cleanup.path, maximum).map_err(E::from)?;
        }
        fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(&cleanup.path)
            .and_then(|file| file.sync_all())
            .map_err(|error| E::from(error.to_string()))?;
        if overwrite {
            fs::rename(&cleanup.path, output)
        } else {
            fs::hard_link(&cleanup.path, output)
        }
        .map_err(|error| E::from(format!("unable to install {}: {error}", output.display())))?;
        Ok(result)
    })()
}

fn protect_inputs(output: &Path, inputs: &[&Path]) -> Result<(), String> {
    let parent = output
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let filename = output.file_name().ok_or("output must name a file")?;
    let destination = fs::canonicalize(parent)
        .map_err(|error| error.to_string())?
        .join(filename);
    let existing = fs::canonicalize(output).ok();
    for input in inputs {
        let input = fs::canonicalize(input).map_err(|error| error.to_string())?;
        if destination == input || existing.as_ref() == Some(&input) {
            return Err("output must not replace an input file or its symlink alias".into());
        }
    }
    Ok(())
}

struct TemporaryOutput {
    path: PathBuf,
    code: &'static str,
}

impl Drop for TemporaryOutput {
    fn drop(&mut self) {
        if let Err(error) = fs::remove_file(&self.path)
            && error.kind() != std::io::ErrorKind::NotFound
        {
            eprintln!(
                "warning[{}]: unable to remove scratch {}: {error}",
                self.code,
                self.path.display()
            );
        }
    }
}

pub(crate) fn check_file_size(path: &Path, maximum: u64) -> Result<(), String> {
    if fs::metadata(path).map_err(|error| error.to_string())?.len() > maximum {
        return Err(format!("HDF5 file exceeds {maximum} bytes"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn publication_rejects_late_destinations_and_recovers_from_owned_write_failures() {
        let root = std::env::temp_dir().join(format!("hdf5-publication-{}", std::process::id()));
        fs::create_dir(&root).unwrap();
        let output = root.join("result.h5");
        let scratch = root.join("result.h5.tmp");
        let error = publish::<(), String>(&output, &[], false, None, "BH0402", |file, _| {
            file.create_group("complete").map_err(hdf5_error)?;
            fs::write(&output, b"arrived during writing").unwrap();
            Ok(())
        })
        .unwrap_err();
        assert!(error.contains("unable to install"));
        assert_eq!(fs::read(&output).unwrap(), b"arrived during writing");
        assert!(!scratch.exists());
        for quota in [None, Some(1)] {
            let error = publish::<(), String>(&output, &[], true, quota, "KR0402", |file, _| {
                file.create_group("incomplete").map_err(hdf5_error)?;
                if quota.is_none() {
                    Err("writer failure".into())
                } else {
                    Ok(())
                }
            })
            .unwrap_err();
            assert!(error.contains(if quota.is_none() {
                "writer failure"
            } else {
                "HDF5 file exceeds"
            }));
            assert!(!scratch.exists());
            assert_eq!(fs::read(&output).unwrap(), b"arrived during writing");
        }
        publish::<(), String>(&output, &[], true, None, "BH0402", |file, _| {
            file.create_group("retry").map_err(hdf5_error)?;
            Ok(())
        })
        .unwrap();
        assert!(File::open(&output).unwrap().group("retry").is_ok());
        assert!(!scratch.exists());
        fs::remove_dir_all(root).unwrap();
    }
}
