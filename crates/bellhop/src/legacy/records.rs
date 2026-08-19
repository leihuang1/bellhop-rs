use std::path::{Path, PathBuf};

use crate::diagnostic::{Diagnostic, SourceLocation};

#[derive(Clone, Debug)]
pub(super) struct Atom {
    pub text: String,
    pub location: SourceLocation,
}

#[derive(Clone, Debug)]
pub(super) struct Slot {
    pub atom: Option<Atom>,
}

#[derive(Clone, Debug)]
enum TokenKind {
    Atom { text: String, quoted: bool },
    Null,
    Slash,
}

#[derive(Clone, Debug)]
struct Token {
    kind: TokenKind,
    location: SourceLocation,
}

pub(super) struct RecordReader<'a> {
    path: PathBuf,
    lines: Vec<&'a str>,
    next_line: usize,
}

impl<'a> RecordReader<'a> {
    pub fn new(source: &'a str, path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            lines: source.lines().collect(),
            next_line: 0,
        }
    }

    pub fn read_fields(
        &mut self,
        field: &'static str,
        expected: usize,
    ) -> Result<Vec<Slot>, Diagnostic> {
        let mut slots = Vec::with_capacity(expected);
        let mut terminated = false;

        while self.next_line < self.lines.len() && slots.len() < expected && !terminated {
            let line_number = self.next_line + 1;
            let line = self.lines[self.next_line];
            self.next_line += 1;
            let tokens = tokenize_line(line, &self.path, line_number, field)?;

            for token in tokens {
                match token.kind {
                    TokenKind::Slash => {
                        terminated = true;
                        break;
                    }
                    TokenKind::Null => {
                        slots.push(Slot { atom: None });
                        if slots.len() >= expected {
                            break;
                        }
                    }
                    TokenKind::Atom { text, quoted } => {
                        append_atom_slots(&mut slots, text, quoted, token.location, expected);
                        if slots.len() >= expected {
                            break;
                        }
                    }
                }
            }
        }

        if slots.len() < expected && !terminated {
            return Err(Diagnostic::error(
                "BH0102",
                format!(
                    "unexpected end of file: expected {expected} value(s), found {}",
                    slots.len()
                ),
                field,
                SourceLocation::new(&self.path, self.lines.len().max(1), 1),
            ));
        }

        Ok(slots)
    }

    pub fn read_required_atom(&mut self, field: &'static str) -> Result<Atom, Diagnostic> {
        let slots = self.read_fields(field, 1)?;
        slots
            .into_iter()
            .next()
            .and_then(|slot| slot.atom)
            .ok_or_else(|| {
                Diagnostic::error(
                    "BH0104",
                    "a value is required here",
                    field,
                    SourceLocation::new(&self.path, self.next_line.max(1), 1),
                )
            })
    }

    pub fn read_string(&mut self, field: &'static str) -> Result<Atom, Diagnostic> {
        self.read_required_atom(field)
    }

    pub fn read_i32(&mut self, field: &'static str) -> Result<(i32, SourceLocation), Diagnostic> {
        let atom = self.read_required_atom(field)?;
        let value = parse_i32(&atom, field)?;
        Ok((value, atom.location))
    }

    pub fn read_f64(&mut self, field: &'static str) -> Result<(f64, SourceLocation), Diagnostic> {
        let atom = self.read_required_atom(field)?;
        let value = parse_f64(&atom, field)?;
        Ok((value, atom.location))
    }
}

pub(super) fn parse_i32(atom: &Atom, field: &'static str) -> Result<i32, Diagnostic> {
    atom.text.parse::<i32>().map_err(|_| {
        Diagnostic::error(
            "BH0103",
            format!("expected an integer, found {:?}", atom.text),
            field,
            atom.location.clone(),
        )
    })
}

pub(super) fn parse_f32(atom: &Atom, field: &'static str) -> Result<f32, Diagnostic> {
    normalized_real(&atom.text).parse::<f32>().map_err(|_| {
        Diagnostic::error(
            "BH0103",
            format!("expected a real number, found {:?}", atom.text),
            field,
            atom.location.clone(),
        )
    })
}

pub(super) fn parse_f64(atom: &Atom, field: &'static str) -> Result<f64, Diagnostic> {
    normalized_real(&atom.text).parse::<f64>().map_err(|_| {
        Diagnostic::error(
            "BH0103",
            format!("expected a real number, found {:?}", atom.text),
            field,
            atom.location.clone(),
        )
    })
}

fn normalized_real(text: &str) -> String {
    text.chars()
        .map(|character| match character {
            'd' => 'e',
            'D' => 'E',
            other => other,
        })
        .collect()
}

fn append_atom_slots(
    slots: &mut Vec<Slot>,
    text: String,
    quoted: bool,
    location: SourceLocation,
    expected: usize,
) {
    if !quoted {
        if let Some((count_text, repeated_text)) = text.split_once('*') {
            if let Ok(count) = count_text.parse::<usize>() {
                for _ in 0..count {
                    if slots.len() >= expected {
                        break;
                    }
                    slots.push(Slot {
                        atom: if repeated_text.is_empty() {
                            None
                        } else {
                            Some(Atom {
                                text: repeated_text.to_owned(),
                                location: location.clone(),
                            })
                        },
                    });
                }
                return;
            }
        }
    }

    slots.push(Slot {
        atom: Some(Atom { text, location }),
    });
}

#[allow(clippy::too_many_lines)]
fn tokenize_line(
    line: &str,
    path: &Path,
    line_number: usize,
    field: &'static str,
) -> Result<Vec<Token>, Diagnostic> {
    let characters: Vec<(usize, char)> = line.char_indices().collect();
    let mut tokens = Vec::new();
    let mut index = 0;
    let mut comma_pending = false;
    let mut has_record_value = false;

    while index < characters.len() {
        let (byte_index, character) = characters[index];
        if character.is_whitespace() {
            index += 1;
            continue;
        }
        if character == ',' {
            if comma_pending || !has_record_value {
                tokens.push(Token {
                    kind: TokenKind::Null,
                    location: SourceLocation::new(
                        path,
                        line_number,
                        line[..byte_index].chars().count() + 1,
                    ),
                });
            }
            comma_pending = true;
            index += 1;
            continue;
        }
        if character == '!' {
            break;
        }

        let location =
            SourceLocation::new(path, line_number, line[..byte_index].chars().count() + 1);
        if character == '/' {
            tokens.push(Token {
                kind: TokenKind::Slash,
                location,
            });
            index += 1;
            continue;
        }
        comma_pending = false;
        has_record_value = true;

        if character == '\'' || character == '"' {
            let quote = character;
            index += 1;
            let mut value = String::new();
            let mut closed = false;
            while index < characters.len() {
                let (_, current) = characters[index];
                if current == quote {
                    if index + 1 < characters.len() && characters[index + 1].1 == quote {
                        value.push(quote);
                        index += 2;
                        continue;
                    }
                    index += 1;
                    closed = true;
                    break;
                }
                value.push(current);
                index += 1;
            }
            if !closed {
                return Err(Diagnostic::error(
                    "BH0101",
                    "unterminated quoted string",
                    field,
                    location,
                ));
            }
            tokens.push(Token {
                kind: TokenKind::Atom {
                    text: value,
                    quoted: true,
                },
                location,
            });
            continue;
        }

        let start = byte_index;
        index += 1;
        while index < characters.len() {
            let current = characters[index].1;
            if current.is_whitespace() || current == ',' || current == '/' || current == '!' {
                break;
            }
            index += 1;
        }
        let end = if index < characters.len() {
            characters[index].0
        } else {
            line.len()
        };
        tokens.push(Token {
            kind: TokenKind::Atom {
                text: line[start..end].to_owned(),
                quoted: false,
            },
            location,
        });
    }

    Ok(tokens)
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::{RecordReader, parse_f64};

    #[test]
    #[allow(clippy::float_cmp)]
    fn understands_quotes_comments_commas_and_slash() {
        let mut reader = RecordReader::new(
            "'a title, with ! text', ! actual comment\n1.0D+3, 2 / ignored\n",
            Path::new("case.env"),
        );

        assert_eq!(
            reader.read_string("title").unwrap().text,
            "a title, with ! text"
        );
        let values = reader.read_fields("values", 4).unwrap();
        assert_eq!(values.len(), 2);
        assert_eq!(
            parse_f64(values[0].atom.as_ref().unwrap(), "values").unwrap(),
            1000.0
        );
        assert_eq!(values[1].atom.as_ref().unwrap().text, "2");
    }

    #[test]
    fn expands_fortran_repetition_and_null_values() {
        let mut reader = RecordReader::new("2*1.5,,3.0, 2* /\n", Path::new("case.env"));
        let values = reader.read_fields("values", 7).unwrap();
        assert_eq!(values.len(), 6);
        assert_eq!(values[0].atom.as_ref().unwrap().text, "1.5");
        assert_eq!(values[1].atom.as_ref().unwrap().text, "1.5");
        assert!(values[2].atom.is_none());
        assert_eq!(values[3].atom.as_ref().unwrap().text, "3.0");
        assert!(values[4].atom.is_none());
        assert!(values[5].atom.is_none());
    }

    #[test]
    fn reports_source_location_for_bad_numbers() {
        let mut reader = RecordReader::new("nope\n", Path::new("case.env"));
        let error = reader.read_f64("frequency").unwrap_err();
        assert_eq!(error.location.line, 1);
        assert_eq!(error.location.column, 1);
        assert_eq!(error.field.as_deref(), Some("frequency"));
    }
}
