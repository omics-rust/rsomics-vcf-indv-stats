use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use rsomics_common::{CommonFlags, Result, ToolMeta};

use rsomics_vcf_indv_stats::{run_depth, run_singletons, run_tstv_summary};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Mode {
    #[value(name = "tstv-summary")]
    TstvSummary,
    #[value(name = "singletons")]
    Singletons,
    #[value(name = "depth")]
    Depth,
}

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-vcf-indv-stats",
    version,
    about = "Per-individual VCF statistics: TsTv-summary, singletons, depth"
)]
pub struct Cli {
    /// VCF/BCF file or - for stdin.
    #[arg(value_name = "VCF")]
    pub vcf: PathBuf,

    /// Statistic to compute.
    #[arg(long = "mode", value_name = "MODE", default_value = "tstv-summary")]
    pub mode: Mode,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let json = self.common.json;
        let path = &self.vcf;

        match self.mode {
            Mode::TstvSummary => {
                let stats = run_tstv_summary(path)?;
                if json {
                    let env = serde_json::json!({
                        "schema_version": rsomics_common::SCHEMA_VERSION,
                        "tool": META.name,
                        "tool_version": META.version,
                        "status": "ok",
                        "result": stats,
                    });
                    println!("{}", serde_json::to_string(&env).unwrap_or_default());
                } else {
                    print!("{}", stats.to_text());
                }
            }
            Mode::Singletons => {
                let scan = run_singletons(path)?;
                // A polyploid site aborts vcftools with exit 1 after it has
                // already written the header and every earlier row. Emit the
                // partial table first, then fail loud so the exit code matches.
                if let Some(site) = scan.abort {
                    if !json {
                        print!("{}", scan.singletons.to_text());
                    }
                    return Err(rsomics_common::RsomicsError::InvalidInput(format!(
                        "Polyploidy found, and not supported by vcftools: {site}"
                    )));
                }
                if json {
                    let env = serde_json::json!({
                        "schema_version": rsomics_common::SCHEMA_VERSION,
                        "tool": META.name,
                        "tool_version": META.version,
                        "status": "ok",
                        "result": scan.singletons,
                    });
                    println!("{}", serde_json::to_string(&env).unwrap_or_default());
                } else {
                    print!("{}", scan.singletons.to_text());
                }
            }
            Mode::Depth => {
                let table = run_depth(path)?;
                if json {
                    let env = serde_json::json!({
                        "schema_version": rsomics_common::SCHEMA_VERSION,
                        "tool": META.name,
                        "tool_version": META.version,
                        "status": "ok",
                        "result": table,
                    });
                    println!("{}", serde_json::to_string(&env).unwrap_or_default());
                } else {
                    print!("{}", table.to_text());
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
