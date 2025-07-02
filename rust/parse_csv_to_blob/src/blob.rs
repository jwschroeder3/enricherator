use std::fs::{File, OpenOptions};
use std::io::{self, Write, Read, BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use crate::BoxedResult;

mod test {
    use super::*;
    use std::collections::HashMap;

    fn set_up_tests(fname: PathBuf) -> (HashMap<String, String>, Blob) {
        let mut map = HashMap::new();
        map.insert("fname".to_string(), "test_data/blob_test/params.blob".to_string());
        map.insert("out_direc".to_string(), "test_data/blob_test".to_string());

        let out_dir = map.get("out_direc").unwrap();
        let blob_file = Blob::create(
            fname,
            400,
            5,
        ).unwrap();

        (map, blob_file)
    }

    #[test]
    fn blob_read() {
        let fname = Path::new("test_data/blob_test").join("params.blob");
        let (map, mut blob_file) = set_up_tests(fname);

        blob_file.read_from_position(0).unwrap();
        assert_eq!(blob_file.buffer.clone(), vec![0,0,0,0,0]);

        blob_file.write_buffer_to_position(1, vec![2]).unwrap();
        blob_file.read_from_position(0).unwrap();
        assert_eq!(blob_file.buffer.clone(), vec![0,2,0,0,0]);

        blob_file.write_buffer_to_position(2, vec![4,5]).unwrap();
        blob_file.read_from_position(1).unwrap();
        assert_eq!(blob_file.buffer.clone(), vec![2,4,5,0,0]);

        assert_eq!(blob_file.cursor_position().unwrap(), 6)
    }

    #[test]
    fn write_f32() {
        let fname = Path::new("test_data/blob_test").join("write_f32.blob");
        let (map, mut blob_file) = set_up_tests(fname);
        let value_to_write: f32 = 12.5;
        // note that, although this is the target buffer, it's not the bytes corresponding
        // to a 12.5 f32 value. That's because we set the buffer size to 4 bytes in the
        // set_up_tests function, so we are going to get one extra byte, which should be 0
        let target = [65, 72, 0, 0, 0];
        blob_file.write_f32_to_position(0, value_to_write).unwrap();
        blob_file.read_from_position(0).unwrap();
        assert_eq!(blob_file.buffer.clone(), target);
    }

}

pub struct Blob {
    file: File,
    pub buffer: Vec<u8>,
}

impl Blob {
    /// creates a binary file of size, "size" bytes
    pub fn create(fname: PathBuf, file_size: usize, buffer_size: usize) -> BoxedResult<Self> {
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .truncate(true)
            .create(true)
            .open(fname)?;
        let file_buffer = vec![0u8; file_size];
        file.write_all(&file_buffer)?;
        let buffer = vec![0u8; buffer_size];
        let blob = Self{file, buffer};
        Ok(blob)
    }

    pub fn open(fname: &str, buffer_size: usize) -> BoxedResult<Self> {
        let file = File::open(fname)?;
        let buffer: Vec<u8> = vec![0; buffer_size];
        Ok(Self{file,buffer})
    }

    pub fn write_value(&mut self, value: f32) -> BoxedResult<()> {
        self.file.write_all(&value.to_be_bytes())?;
        Ok(())
    }

    pub fn read_from_position(
            &mut self,
            start: u64,
    ) -> BoxedResult<()> {
        // Seek to the desired offset position from the start byte
        self.move_cursor(start)?;
        // Read the exact number of bytes into the buffer
        self.file.read_exact(&mut self.buffer)?;
        // after this action, file cursor has moved `buffer size` bytes
        Ok(())
    }

    pub fn write_buffer_to_position(
            &mut self,
            start: u64,
            buffer: Vec<u8>,
    ) -> BoxedResult<()> {
        // Seek to the desired offset position from the start byte
        self.move_cursor(start)?;
        self.file.write_all(&buffer)?;
        Ok(())
    }

    pub fn write_f32_to_position(
            &mut self,
            start: u64,
            value: f32,
    ) -> BoxedResult<()> {
        // Seek to the desired offset position from the start byte
        self.move_cursor(start)?;
        let bytes = value.to_be_bytes();
        self.file.write_all(&bytes)?;
        Ok(())
    }

    pub fn cursor_position(&mut self) -> BoxedResult<u64> {
        let pos = self.file.seek(SeekFrom::Current(0))?;
        Ok(pos)
    }

    fn move_cursor(&mut self, offset: u64) -> BoxedResult<()> {
        self.file.seek(SeekFrom::Start(offset))?;
        Ok(())
    }

    pub fn read_f32_from_position(&mut self, offset: u64) -> BoxedResult<Vec<f32>> {
        self.read_from_position(offset)?;
        let mut vals: Vec<f32> = Vec::new();
        for start_idx in (0..self.buffer.len()).step_by(4) {
            let end_idx = start_idx + 4;
            let val_bytes: [u8; 4] =  self.buffer[start_idx..end_idx].try_into()?;
            let val = f32::from_be_bytes(val_bytes);
            vals.push(val);
        }
        Ok(vals)
    }
}

