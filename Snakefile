from snakemake.utils import min_version
min_version("5.3.0")

configfile: "config.yml"

rule raw_data_process_generate:
    input:
        report = "code/raw_data_process.Rmd" 
    output:
        raw_report = "raw_data_process.html", 
        cell = "data/processed/cells_{chemical}_{cell_line}.csv",
        info = "data/processed/info_{chemical}_{cell_line}.csv"
    params:
        cell_line = config["cell_line"],
        chemical = config["chemical"],
        outdir = "data/processed"
    shell:
        """
        Rscript -e 'parameters <- list(root_directory = getwd(),
                                       cell_line = "{params.cell_line}",
                                       chemical = "{params.chemical}",
                                       cell = "{output.cell}",
                                       info = "{output.info}");
                    rmarkdown::render("{input.report}", 
                                      output_file= "{output.raw_report}", 
                                      output_dir="{params.outdir}", 
                                      params= parameters)' 
        """
rule raw_data_process_generate:
    input:
        report = "code/raw_data_process.Rmd" 
    output:
        "raw_data_process.html", 
        expand("data/processed/cells_{chemical}_{cell_line}.csv",cell_line = config["cell_line"],chemical = config["chemical"])
    params:
        cell_line = config["cell_line"],
        chemical = config["chemical"],
        outdir = "data/processed"
    shell:
        """
        Rscript -e 'parameters <- list(root_directory = getwd(),
                                       cell_line = "{params.cell_line}",
                                       chemical = "{params.chemical}");
                    rmarkdown::render("{input.report}", 
                                      output_file= "{output}", 
                                      output_dir="{params.outdir}", 
                                      params= parameters)' 
        """

rule qc:
    input:
        report = "code/QC.Rmd",
        cells = "data/processed/cells_{chemical}_{cell_line}.csv",
        info = "data/processed/info_{chemical}_{cell_line}.csv"
    output:
        "QC_{chemical}_{cell_line}.html"
    params:
        outdir = "results/reports"
    shell:
        """
        Rscript -e 'parameters <- list(root_directory = getwd(),
                                       cells = "{input.cells}",
                                       info = "{input.info}");
                    rmarkdown::render("{input.report}", 
                                      output_file= "{output}", 
                                      output_dir="{params.outdir}", 
                                      params= parameters)' 
        """

rule DE_analysis:
    input:
        report = "code/DE_analysis.Rmd",
        cells = "data/processed/cells_{chemical}_{cell_line}.csv",
        info = "data/processed/info_{chemical}_{cell_line}.csv"
    output:
        output_report = "DE_analysis_{group}_{cell_line}.html", 
        DE_gene = "sig_0_05_DE_genes_{group}_{cell_line}.csv",
        all_gene = "background_genes_shrunk_{group}_{cell_line}.csv"
    params:
        contrast = config["contrast"],
        group = config["group"],
        outdir = "results/reports"
    shell:
        """
        Rscript -e 'parameters <- list(root_directory = getwd(),
                                       cells = "{input.cells}",
                                       info = "{input.info}",
                                       group = "{input.group}",
                                       contrast = "{input.contrast}",
                                       DE_gene = "{output.DE_gene}",
                                       all_gene = "{output.all_gene}");
                    rmarkdown::render("{input.report}", 
                                      output_file= "{output.output_report}", 
                                      output_dir="{params.outdir}", 
                                      params= parameters)' 
        """

rule functional_analysis:
    input:
        report = "code/functional_analysis.Rmd",
        DE_gene = "results/tables/DES/sig_0_05_DE_genes_{group}_{cell_line}.csv",
        all_gene = "results/tables/DES/background_genes_shrunk_{group}_{cell_line}.csv"
    output:
        GO_report = "DE_analysis_{group}_{cell_line}.html"
    params:
        contrast = config["contrast"],
        group = config["group"],
        outdir = "results/reports"
    shell:
        """
        Rscript -e 'parameters <- list(root_directory = getwd(),
                                       group = "{input.group}",
                                       contrast = "{input.contrast}");
                    rmarkdown::render("{input.report}", 
                                      output_file= "{output.GO_report}", 
                                      output_dir="{params.outdir}", 
                                      params= parameters)' 
        """