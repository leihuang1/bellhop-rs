//! Shared HDF5 metadata I/O, not a shared schema.

use hdf5::H5Type;
use hdf5::types::VarLenUnicode;
use std::str::FromStr;

pub(crate) fn write_scalar_attribute<T: H5Type>(
    parent: &hdf5::Location,
    name: &str,
    value: &T,
) -> Result<(), String> {
    parent
        .new_attr::<T>()
        .create(name)
        .and_then(|attribute| attribute.write_scalar(value))
        .map_err(hdf5_error)
}

pub(crate) fn write_string_attribute(
    parent: &hdf5::Location,
    name: &str,
    value: &str,
) -> Result<(), String> {
    let value = VarLenUnicode::from_str(value)
        .map_err(|error| format!("invalid metadata string: {error}"))?;
    parent
        .new_attr::<VarLenUnicode>()
        .create(name)
        .and_then(|attribute| attribute.write_scalar(&value))
        .map_err(hdf5_error)
}

pub(crate) fn write_string_array_attribute(
    parent: &hdf5::Location,
    name: &str,
    values: &[String],
) -> Result<(), String> {
    let values = values
        .iter()
        .map(|value| {
            VarLenUnicode::from_str(value)
                .map_err(|error| format!("invalid metadata string: {error}"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    parent
        .new_attr_builder()
        .with_data(&values)
        .create(name)
        .map(|_| ())
        .map_err(hdf5_error)
}

#[allow(clippy::needless_pass_by_value)]
pub(crate) fn hdf5_error(error: hdf5::Error) -> String {
    format!("HDF5 output failed: {error}")
}
