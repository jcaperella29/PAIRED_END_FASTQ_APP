library(shiny)
library(shinythemes)
library(DT)
library(ShortRead)
library(shinyFiles)
library(fs)
library(future)
library(future.apply)
library(Rsubread)
library(promises)
library(ballgown)
library(rtracklayer)
library(seqTools)
library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(GenomicAlignments)
library(Biostrings)

# Define UI
ui <- fluidPage(
  theme = shinytheme("cyborg"),
  titlePanel("JCap Paired_End FASTQ Processing App"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("directory", "Set Working Directory", "Please select a directory"),
      fileInput("fastq_file1", "Upload FASTQ File 1", accept = c(".fastq", ".fq")),
      fileInput("fastq_file2", "Upload FASTQ File 2", accept = c(".fastq", ".fq")),
      selectInput("machine_type", "Select Sequencing Machine Type", 
                  choices = c("Illumina/Sanger (Phred+33)" = "phred33", 
                              "Solexa (Phred+64)" = "phred64")),
      textInput("sample_name", "Name of Sample", value = "Sample1"),
      fileInput("barcodes_file", "Upload Barcodes File (Optional)", accept = c(".csv", ".tsv", ".txt")),
      actionButton("demultiplex", "Run Demultiplexing"),
      fileInput("demultiplexed_file1", "Upload Demultiplexed File 1", accept = c(".fastq", ".fq")),
      fileInput("demultiplexed_file2", "Upload Demultiplexed File 2", accept = c(".fastq", ".fq")),
      actionButton("generate_summary", "Generate FASTQ Summary"),
      downloadButton("download_summary", "Download Summary as CSV"),
      sliderInput("quality_threshold", "Quality Threshold", min = 0, max = 100, value = 20),
      actionButton("filter_reads", "Filter Reads Below Quality Threshold"),
      fileInput("filtered_fastq_file1", "Upload Filtered FASTQ File 1 (Optional)", accept = c(".fastq", ".fq")),
      fileInput("filtered_fastq_file2", "Upload Filtered FASTQ File 2 (Optional)", accept = c(".fastq", ".fq")),
      fileInput("genome_fasta_indexing", "Upload Genome FASTA File for Indexing", accept = c(".fasta", ".fa")),
      numericInput("memory_size", "Memory Allocation (MB)", value = 8000, min = 1000, step = 1000),
      textInput("index_base_name", "Index Base Name for Alignment (e.g., hg38_index)", value = "custom_genome_index"),
      actionButton("create_index", "Create Genome Index"),
      actionButton("align_reads", "Align FASTQs to Genome"),
      fileInput("bam_file", "Upload BAM File", accept = c(".bam")),
      actionButton("convert_bam_to_bed", "Convert BAM to BED"),
      fileInput("gtf_file", "Upload GTF File", accept = c(".gtf")),
      actionButton("identify_gtf_attributes", "Identify GTF Attributes"),
      textInput("gtf_attr_type", "GTF Attribute Type (e.g., gene_id)", value = "gene_id"),
      actionButton("generate_counts_matrix", "Generate Counts Matrix"),
      downloadButton("download_counts_matrix", "Download Counts Matrix"),
      actionButton("create_fai", "Generate FAI"),
      actionButton("create_gff3", "Generate GFF3")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("FASTQ Summary", DTOutput("fastq_summary")),
        tabPanel("Counts Matrix", DTOutput("counts_matrix")),
        tabPanel("GTF Attributes", verbatimTextOutput("gtf_attributes")),
        tabPanel("README", verbatimTextOutput("readme_contents"))
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 2 * 1024^50000)
  
  readme_file_path <- "JCAP Paired_End FASTQ Processing App_ReadMe.txt" # Update this path to your actual README file
  
  # Reactive expression to read the README file
  readme_contents <- reactive({
    if (file.exists(readme_file_path)) {
      readLines(readme_file_path, warn = FALSE)
    } else {
      "README file not found."
    }
  })
  
  # Render the README content
  output$readme_contents <- renderText({
    paste(readme_contents(), collapse = "\n")
  })
  
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  shinyDirChoose(input, "directory", roots = volumes, session = session)
  
  working_dir <- reactiveVal(getwd())
  
  observeEvent(input$directory, {
    req(input$directory)
    dir_selected <- parseDirPath(volumes, input$directory)
    dir_selected <- as.character(dir_selected)
    
    if (length(dir_selected) > 0) {
      setwd(dir_selected)
      working_dir(dir_selected)
      showNotification(paste("Working directory set to:", dir_selected))
    } else {
      showNotification("No directory selected", type = "error")
    }
  })
  
  fastq_files <- reactiveVal(list())
  filtered_files <- reactiveVal(list())
  demultiplexed_files <- reactiveVal(list(file1 = NULL, file2 = NULL))
  counts_matrix <- reactiveVal(NULL)
  bed_file_path <- reactiveVal(NULL)
  genome_fasta_path_indexing <- reactiveVal(NULL)
  
  observeEvent(input$fastq_file1, {
    req(input$fastq_file1)
    showNotification("FASTQ File 1 uploaded", duration = 2)
    fastq_files(list(file1 = input$fastq_file1, file2 = fastq_files()$file2))
    showNotification("FASTQ File 1 uploaded successfully.")
  })
  
  observeEvent(input$fastq_file2, {
    req(input$fastq_file2)
    showNotification("FASTQ File 2 uploaded", duration = 2)
    fastq_files(list(file1 = fastq_files()$file1, file2 = input$fastq_file2))
    showNotification("FASTQ File 2 uploaded successfully.")
  })
  
  observeEvent(input$filtered_fastq_file1, {
    req(input$filtered_fastq_file1)
    showNotification("Filtered FASTQ File 1 uploaded", duration = 2)
    filtered_files(list(file1 = input$filtered_fastq_file1, file2 = filtered_files()$file2))
    showNotification("Filtered FASTQ File 1 uploaded successfully.")
  })
  
  observeEvent(input$filtered_fastq_file2, {
    req(input$filtered_fastq_file2)
    showNotification("Filtered FASTQ File 2 uploaded", duration = 2)
    filtered_files(list(file1 = filtered_files()$file1, file2 = input$filtered_fastq_file2))
    showNotification("Filtered FASTQ File 2 uploaded successfully.")
  })
  
  observeEvent(input$demultiplexed_file1, {
    req(input$demultiplexed_file1)
    showNotification("Demultiplexed File 1 uploaded", duration = 2)
    demultiplexed_files(list(file1 = input$demultiplexed_file1, file2 = demultiplexed_files()$file2))
    showNotification("Demultiplexed File 1 uploaded successfully.")
  })
  
  observeEvent(input$demultiplexed_file2, {
    req(input$demultiplexed_file2)
    showNotification("Demultiplexed File 2 uploaded", duration = 2)
    demultiplexed_files(list(file1 = demultiplexed_files()$file1, file2 = input$demultiplexed_file2))
    showNotification("Demultiplexed File 2 uploaded successfully.")
  })
  
  observeEvent(input$genome_fasta_indexing, {
    req(input$genome_fasta_indexing)
    showNotification("Genome FASTA File for Indexing uploaded", duration = 2)
    genome_fasta_path_indexing(input$genome_fasta_indexing$datapath)
  })
  
  summary_df <- reactiveVal(NULL)  # Store the summary data frame for download
  
  observeEvent(input$generate_summary, {
    req(
      (isTruthy(input$fastq_file1) && isTruthy(input$fastq_file2)) ||
        (isTruthy(demultiplexed_files()$file1) && isTruthy(demultiplexed_files()$file2)) ||
        (isTruthy(filtered_files()$file1) && isTruthy(filtered_files()$file2))
    )
    
    filename1 <- if (!isTruthy(filtered_files()$file1)) {
      if (!isTruthy(demultiplexed_files()$file1)) {
        input$fastq_file1$name
      } else {
        demultiplexed_files()$file1$name
      }
    } else {
      filtered_files()$file1$name
    }
    
    filename2 <- if (!isTruthy(filtered_files()$file2)) {
      if (!isTruthy(demultiplexed_files()$file2)) {
        input$fastq_file2$name
      } else {
        demultiplexed_files()$file2$name
      }
    } else {
      filtered_files()$file2$name
    }
    
    file_to_use1 <- if (isTruthy(filtered_files()$file1)) {
      filtered_files()$file1$datapath
    } else if (isTruthy(demultiplexed_files()$file1)) {
      demultiplexed_files()$file1$datapath
    } else {
      input$fastq_file1$datapath
    }
    
    file_to_use2 <- if (isTruthy(filtered_files()$file2)) {
      filtered_files()$file2$datapath
    } else if (isTruthy(demultiplexed_files()$file2)) {
      demultiplexed_files()$file2$datapath
    } else {
      input$fastq_file2$datapath
    }
    
    req(file_to_use1, file_to_use2)  # Both files must be present
    
    showNotification("Generating FASTQ summary", duration = 2)
    
    quality_threshold <- input$quality_threshold
    phred_offset <- get_phred_offset(input$machine_type)
    
    future({
      quality_scores1 <- extract_quality_scores(file_to_use1)
      quality_scores2 <- extract_quality_scores(file_to_use2)
      
      res1 <- process_chunk(quality_scores1, quality_threshold)
      res2 <- process_chunk(quality_scores2, quality_threshold)
      
      summary_df <- data.frame(
        Filename = c(filename1, filename2),
        NumReads = c(res1$num_reads, res2$num_reads),
        PercentAboveThreshold = c(res1$percent_above_threshold, res2$percent_above_threshold),
        AvgQuality = c(res1$avg_quality, res2$avg_quality)
      )
      summary_df
    }, seed = TRUE) %...>% (function(summary_df_result) {
      summary_df(summary_df_result)  # Store summary data for download
      output$fastq_summary <- renderDT({
        datatable(summary_df_result)
      })
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error generating FASTQ summary:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  output$download_summary <- downloadHandler(
    filename = function() {
      paste("FASTQ_summary-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(summary_df())  # Ensure summary exists before allowing download
      write.csv(summary_df(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$create_index, {
    req(input$genome_fasta_indexing, input$memory_size)
    
    genome_fasta <- input$genome_fasta_indexing$datapath
    memory_size <- input$memory_size
    working_directory <- working_dir()
    
    showNotification("Starting genome index creation. This may take between 1.5 to 2 hours.", duration = NULL)
    
    future({
      setwd(working_directory)
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      index_base <- paste0("index_", timestamp)
      
      buildindex(
        basename = index_base,
        reference = genome_fasta,
        gappedIndex = FALSE,
        indexSplit = TRUE,
        memory = memory_size,
        TH_subread = 100,
        colorspace = FALSE
      )
      
      index_files <- list.files(pattern = paste0(index_base, ".*\\.(reads|files|txt|b.array|b.tab)"), path = working_directory, full.names = TRUE)
      
      list(index_files = index_files, index_base = index_base)
    }, seed = TRUE) %...>% (function(result) {
      index_files <- result$index_files
      if (length(index_files) > 0) {
        showNotification("Index created successfully.", duration = 2)
        output$index_creation_log <- renderText(paste("Index created successfully:", paste(index_files, collapse = "\n")))
        output$index_file_path <- renderText({
          paste("Generated Index File Paths: ", paste(index_files, collapse = ", "))
        })
      } else {
        showNotification("Index creation failed.", type = "error", duration = 2)
        output$index_creation_log <- renderText("Index creation failed.")
      }
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error creating index:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  observeEvent(input$align_reads, {
    files_to_use <- if (!is.null(filtered_files()$file1) && !is.null(filtered_files()$file2)) {
      c(filtered_files()$file1$datapath, filtered_files()$file2$datapath)
    } else if (!is.null(demultiplexed_files()$file1) && !is.null(demultiplexed_files()$file2)) {
      c(demultiplexed_files()$file1$datapath, demultiplexed_files()$file2$datapath)
    } else {
      c(fastq_files()$file1$datapath, fastq_files()$file2$datapath)
    }
    index_base_name <- input$index_base_name
    
    req(index_base_name)
    showNotification("Starting alignment of FASTQs to genome. This may take between 1 to 2 hours.", duration = NULL)
    
    working_directory <- working_dir()
    
    future({
      setwd(working_directory)
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      bam_file <- file.path(working_directory, paste0("aligned_", timestamp, ".bam"))
      
      align(index = index_base_name, 
            readfile1 = files_to_use[1], 
            readfile2 = if (length(files_to_use) > 1) files_to_use[2] else NULL, 
            output_file = bam_file)
      bam_file
    }, seed = TRUE) %...>% (function(bam_file) {
      showNotification("Alignment completed.", duration = 2)
      output$align_log <- renderText(paste("Alignment completed. BAM file: ", bam_file))
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error during alignment:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  output$download_counts_matrix <- downloadHandler(
    filename = function() {
      paste("counts_matrix-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(counts_matrix(), file, row.names = TRUE)
    }
  )
  
  observeEvent(input$generate_counts_matrix, {
    req(input$gtf_file, input$bam_file, input$gtf_attr_type, input$sample_name)
    showNotification("Generating counts matrix", duration = 2)
    
    gtf_file <- input$gtf_file$datapath
    bam_file <- input$bam_file$datapath
    working_directory <- working_dir()
    gtf_attr_type <- input$gtf_attr_type
    sample_name <- isolate(input$sample_name)
    
    future({
      temp_dir <- file.path(working_directory, "temp_featureCounts")
      dir.create(temp_dir, showWarnings = FALSE)
      
      setwd(temp_dir)
      
      on.exit({
        setwd(working_directory)
        unlink(temp_dir, recursive = TRUE)
      }, add = TRUE)
      
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      tryCatch({
        result <- featureCounts(
          files = bam_file,
          annot.ext = gtf_file,
          isGTFAnnotationFile = TRUE,
          GTF.attrType = gtf_attr_type,
          isPairedEnd = TRUE
        )
        counts_df <- as.data.frame(result$counts)
        colnames(counts_df) <- sample_name
        output_file <- file.path(working_directory, paste0("counts_", timestamp, ".csv"))
        write.csv(counts_df, output_file, row.names = TRUE)
        list(counts_matrix = counts_df, output_file = output_file, success = TRUE)
      }, error = function(e) {
        list(success = FALSE, message = e$message)
      })
    }, seed = TRUE) %...>% (function(result) {
      if (result$success) {
        counts_matrix(result$counts_matrix)
        output$counts_matrix <- renderDT({
          datatable(result$counts_matrix)
        })
        showNotification(paste("Counts matrix generation completed. File saved to:", result$output_file), duration = 2)
      } else {
        showModal(modalDialog(
          title = "Error",
          paste("Error generating counts matrix:", result$message),
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error generating counts matrix:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  observeEvent(input$demultiplex, {
    req(input$barcodes_file, input$fastq_file1)
    
    barcodes_file <- input$barcodes_file$datapath
    fastq_file1 <- input$fastq_file1$datapath
    output_dir <- working_dir()
    
    future({
      barcodes <- read.csv(barcodes_file, header = FALSE, stringsAsFactors = FALSE, skip = 1)
      barcodes <- as.character(barcodes[[1]])
      barcodes <- gsub("\\s", "", barcodes)
      
      barcode_lengths <- unique(nchar(barcodes))
      if (length(barcode_lengths) > 1) {
        stop("All barcodes must be of the same length.")
      }
      
      fq <- readFastq(fastq_file1)
      
      sequences <- as.character(sread(fq))
      qualities <- extract_quality_scores(fastq_file1)
      
      demultiplexed <- custom_demultiplex(fastq_file1, barcodes, output_dir)
      
      for (barcode in names(demultiplexed)) {
        demux_reads <- demultiplexed[[barcode]]
        
        if (length(demux_reads$sequences) > 0) {
          demux_fq <- ShortReadQ(
            sread = DNAStringSet(demux_reads$sequences),
            quality = new(Class = class(quality(fq)), demux_reads$qualities),
            id = BStringSet(paste0("@SEQ_ID_", 1:length(demux_reads$sequences)))
          )
          writeFastq(demux_fq, file.path(output_dir, paste0("demultiplexed_", barcode, ".fastq")), mode = "w", compress = FALSE, append = FALSE)
        }
      }
      TRUE
    }, seed = TRUE) %...>% (function(result) {
      showNotification("Demultiplexing completed.")
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error during demultiplexing:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  observeEvent(input$create_fai, {
    req(input$genome_fasta_indexing)
    showNotification("Generating FAI file", duration = 2)
    
    fasta_file <- input$genome_fasta_indexing$datapath
    working_directory <- working_dir()
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    future({
      fai_file <- file.path(working_directory, paste0("genome_", timestamp, ".fai"))
      seqTools::writeFai(fasta_file, fai_file)
      fai_file
    }, seed = TRUE) %...>% (function(fai_file) {
      showNotification("FAI file generated successfully.", duration = 2)
      output$fai_file_path <- renderText({
        paste("Generated FAI File Path: ", fai_file)
      })
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error generating FAI file:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  observeEvent(input$create_gff3, {
    req(input$gtf_file)
    showNotification("Generating GFF3 file", duration = 2)
    
    gtf_file <- input$gtf_file$datapath
    working_directory <- working_dir()
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    future({
      gff_file <- file.path(working_directory, paste0("annotations_", timestamp, ".gff3"))
      gtf <- import(gtf_file, format = "gtf")
      export(gtf, gff_file, format = "gff3")
      gff_file
    }, seed = TRUE) %...>% (function(gff_file) {
      showNotification("GFF3 file generated successfully.", duration = 2)
      output$gff3_file_path <- renderText({
        paste("Generated GFF3 File Path: ", gff_file)
      })
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error generating GFF3 file:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  observeEvent(input$identify_gtf_attributes, {
    req(input$gtf_file)
    showNotification("Identifying GTF attributes", duration = 2)
    
    gtf_file <- input$gtf_file$datapath
    
    future({
      identify_gtf_attributes(gtf_file)
    }, seed = TRUE) %...>% (function(attributes) {
      output$gtf_attributes <- renderText(paste("Identified GTF attributes:", paste(attributes, collapse = ", ")))
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error identifying GTF attributes:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  observeEvent(input$convert_bam_to_bed, {
    req(input$bam_file)
    showNotification("Converting BAM to BED", duration = 2)
    
    bam_file <- input$bam_file$datapath
    working_directory <- working_dir()
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    bed_file <- file.path(working_directory, paste0("converted_", timestamp, ".bed"))
    
    future({
      convert_bam_to_bed(bam_file, bed_file)
    }, seed = TRUE) %...>% (function(bed_file) {
      bed_file_path(bed_file)
      showNotification("BAM to BED conversion completed.", duration = 2)
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error during BAM to BED conversion:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
