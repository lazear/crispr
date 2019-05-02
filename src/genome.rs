//! Utilities for loading genomic information
use memchr::{memchr_iter, Memchr};
use std::collections::{BTreeMap, HashMap};
use std::fs::{File, OpenOptions};
use std::io::{self, prelude::*};
use std::path::Path;
use std::str;

struct Pitchfork<'a> {
    pos: usize,
    haystack: &'a [u8],
    inner: Memchr<'a>,
}

impl<'a> Pitchfork<'a> {
    pub fn new(needle: u8, haystack: &'a [u8]) -> Self {
        Self {
            pos: 0,
            haystack,
            inner: memchr_iter(needle, haystack),
        }
    }
}

impl<'a> Iterator for Pitchfork<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let end = match self.inner.next() {
            Some(e) => e,
            None => {
                if self.pos < self.haystack.len() {
                    self.haystack.len()
                } else {
                    return None;
                }
            }
        };
        let slice = &self.haystack[self.pos..end];
        self.pos = end + 1;
        Some(slice)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone)]
pub struct Genome {
    pub map: HashMap<String, Span>,
    condensed: String,
    index: BTreeMap<usize, String>,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Entry<'a> {
    pub id: &'a str,
    pub sequence: &'a str,
    pub span: Span,
}

impl Genome {
    /// Build a reference genome from a Fasta file
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Genome> {
        let mut buf = String::new();
        File::open(path)?.read_to_string(&mut buf)?;

        let mut map = HashMap::new();
        let mut index = BTreeMap::new();
        let mut condensed = String::with_capacity(buf.len());

        let mut length = 0usize;
        let mut start = 0usize;

        let mut iter = Pitchfork::new('\n' as u8, buf.as_bytes());
        let mut last = iter.next().unwrap();

        for line in iter {
            if line[0] == '>' as u8 {
                let id = str::from_utf8(last).unwrap();
                last = line;

                let enst = &id[1..16];
                let span = Span { start, end: length };
                start = length;
                assert_eq!(length, condensed.len());

                index.insert(span.start, enst.into());
                map.insert(enst.into(), span);
            } else {
                length += line.len();
                condensed.push_str(str::from_utf8(line).unwrap());
            }
        }

        Ok(Genome {
            map,
            index,
            condensed,
        })
    }

    /// Given an index into the condensed genome, return the entry corresponding
    /// to the closest previous entry in the Fasta file
    pub fn range(&self, pos: usize) -> Option<Entry<'_>> {
        let (_, id) = self.index.range(..pos).next_back()?;
        let span = *self.map.get(id)?;

        Some(Entry {
            id,
            sequence: &self.condensed[span.start..span.end],
            span,
        })
    }

    /// Return a reference to the condensed genome
    pub fn condensed(&self) -> &str {
        self.condensed.as_ref()
    }

    pub fn get_by_ensembl(&self, id: &str) -> Option<&str> {
        let span = *self.map.get(id)?;
        Some(&self.condensed[span.start..span.end])
    }

    /// Write to a file located at `path` - this is essentially a flattened
    /// Fasta format, where each sequence is on one file only.
    pub fn write<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut f = OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(true)
            .open(path)?;
        for (id, span) in &self.map {
            let seq = &self.condensed[span.start..span.end];
            writeln!(f, "{}\n{}", id, seq)?;
        }
        Ok(())
    }
}
