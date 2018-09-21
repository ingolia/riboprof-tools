extern crate failure;

extern crate riboprof;

use std::io;
use std::io::Write;
use std::process;

use riboprof::fastx_split::*;

fn main() {
    match wrapper() {
        Err(e) => {
            io::stderr().write(format!("{}\n", e).as_bytes()).unwrap();
            process::exit(1);
        }
        _ => (),
    };
}

fn wrapper() -> Result<(), failure::Error> {
    let cli = CLI::new()?;
    let config = Config::new(&cli)?;
    fastx_split(config)
}
