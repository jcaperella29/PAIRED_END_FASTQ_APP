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

# Set up parallel processing
plan(multisession, workers = 2)

# Functions to generate FAI and GFF3 files
generate_fai <- function(fasta_file, output_dir) {
  require(seqTools)
  fai_file <- file.path(output_dir, paste0(basename(fasta_file), ".fai"))
  cat("Generating FAI file at:", fai_file, "\n")
  seqTools::writeFai(fasta_file, fai_file)
  fai_file
}

generate_gff3 <- function(gtf_file, output_dir) {
  gff_file <- file.path(output_dir, sub("\\.gtf$", ".gff3", basename(gtf_file)))
  cat("Generating GFF3 file at:", gff_file, "\n")
  gtf <- import(gtf_file, format = "gtf")
  export(gtf, gff_file, format = "gff3")
  gff_file
}

# Function to convert BAM to BED
convert_bam_to_bed <- function(bam_file, bed_file) {
  bam <- import(bam_file, format = "BAM")
  if (length(bam) == 0) {
    stop("The BAM file is empty.")
  }
  
  bed <- data.frame(
    chrom = as.character(seqnames(bam)),
    chromStart = start(bam) - 1,
    chromEnd = end(bam)
  )
  
  bed <- unique(bed)
  
  write.table(bed, file = bed_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  return(bed_file)
}

# Function to print BAM header and content
print_bam_content <- function(bam_file) {
  header <- scanBamHeader(bam_file)
  print(header)
  param <- ScanBamParam(what = scanBamWhat())
  bam <- scanBam(bam_file, param = param)
  print(bam)
}

# Function to get Phred offset based on machine type
get_phred_offset <- function(machine_type) {
  switch(machine_type,
         "phred33" = 33,
         "phred64" = 64,
         stop("Unknown machine type"))
}

filter_fastq_by_quality <- function(fastq_file, output_dir, phred_offset, quality_threshold) {
  # Generate a timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Read the FASTQ sequences and their quality scores
  fq <- readFastq(fastq_file)
  quality_scores <- extract_quality_scores(fastq_file, phred_offset)
  
  # Identify reads that meet the quality threshold
  reads_above_threshold <- which(sapply(quality_scores, function(scores) {
    mean(scores) >= quality_threshold
  }))
  
  # Filter sequences and quality scores based on the threshold
  filtered_seqs <- fq[reads_above_threshold]
  
  # Generate the output filename with timestamp
  output_file <- file.path(output_dir, paste0("filtered_", basename(fastq_file), "_", timestamp, ".fastq"))
  
  # Write the filtered sequences to the new FASTQ file, overwriting if it exists
  writeFastq(filtered_seqs, output_file, mode = "w", compress = FALSE, append = FALSE)
  
  return(output_file)
}


# Function to read a FASTQ file and extract quality scores with Phred offset consideration
extract_quality_scores <- function(fastq_file, phred_offset) {
  con <- file(fastq_file, open = "r")
  quality_scores <- list()
  line_counter <- 0
  current_quality_scores <- NULL
  
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line_counter <- line_counter + 1
    if (line_counter %% 4 == 0) {
      current_quality_scores <- strsplit(line, "")[[1]]
      quality_scores <- c(quality_scores, list(current_quality_scores))
    }
  }
  
  close(con)
  
  convert_to_phred <- function(ascii_chars) {
    sapply(ascii_chars, function(char) {
      as.integer(charToRaw(char)) - phred_offset
    })
  }
  
  quality_scores <- lapply(quality_scores, convert_to_phred)
  return(quality_scores)
}

# Function to convert numeric Phred scores to ASCII characters
convert_numeric_to_ascii <- function(numeric_qualities, phred_offset) {
  sapply(numeric_qualities, function(score) {
    rawToChar(as.raw(score + phred_offset))
  })
}


# Function to process quality scores and compute statistics
process_chunk <- function(quality_scores, quality_threshold) {
  num_reads <- length(quality_scores)
  avg_quality <- mean(unlist(quality_scores))
  percent_above_threshold <- mean(unlist(quality_scores) >= quality_threshold) * 100
  
  list(
    avg_quality = avg_quality,
    percent_above_threshold = percent_above_threshold,
    num_reads = num_reads
  )
}

filter_fastq_by_quality <- function(fastq_file, output_dir, phred_offset, quality_threshold) {
  # Generate a timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Read the FASTQ sequences and their quality scores
  fq <- readFastq(fastq_file)
  quality_scores <- extract_quality_scores(fastq_file, phred_offset)
  
  # Identify reads that meet the quality threshold
  reads_above_threshold <- which(sapply(quality_scores, function(scores) {
    mean(scores) >= quality_threshold
  }))
  
  # Filter sequences and quality scores based on the threshold
  filtered_seqs <- fq[reads_above_threshold]
  
  # Generate the output filename with timestamp
  output_file <- file.path(output_dir, paste0("filtered_", basename(fastq_file), "_", timestamp, ".fastq"))
  
  # Write the filtered sequences to the new FASTQ file
  writeFastq(filtered_seqs, output_file, mode = "w", compress = FALSE, append = FALSE)
  
  return(output_file)
}

custom_demultiplex <- function(fastq_file, barcodes, output_dir, phred_offset) {
  sequences <- sread(readFastq(fastq_file))
  qualities <- extract_quality_scores(fastq_file, phred_offset)  # Pass phred_offset here
  
  demux_results <- list()
  
  for (barcode in barcodes) {
    matched_indices <- grep(barcode, as.character(sequences))
    if (length(matched_indices) > 0) {
      demux_results[[barcode]] <- list(
        sequences = as.character(sequences[matched_indices]),
        qualities = qualities[matched_indices]
      )
    }
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  for (barcode in names(demux_results)) {
    output_file <- file.path(output_dir, paste0("demultiplexed_", barcode, "_", timestamp, ".fastq"))
    seqs <- demux_results[[barcode]]$sequences
    quals <- demux_results[[barcode]]$qualities
    
    quality_scores <- BStringSet(sapply(quals, function(q) {
      paste(convert_numeric_to_ascii(q, phred_offset), collapse = "")
    }))
    
    # Overwrite existing file if it exists
    writeFastq(ShortReadQ(sread = DNAStringSet(seqs), quality = quality_scores), output_file, mode = "w", compress = FALSE, append = FALSE)
  }
  
  return(output_dir)
}

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
  observeEvent(input$filter_reads, {
    req(fastq_files()$file1, fastq_files()$file2)
    
    # Notify the user that filtering has started
    showNotification("Filtering reads below quality threshold started.", duration = 2)
    
    future({
      phred_offset <- get_phred_offset(input$machine_type)
      quality_threshold <- input$quality_threshold
      
      # Define the output directory
      output_dir <- working_dir()
      
      # Filter the reads for each FASTQ file with a timestamped name
      filtered_file1 <- filter_fastq_by_quality(fastq_files()$file1$datapath, output_dir, phred_offset, quality_threshold)
      filtered_file2 <- filter_fastq_by_quality(fastq_files()$file2$datapath, output_dir, phred_offset, quality_threshold)
      
      list(filtered_file1 = filtered_file1, filtered_file2 = filtered_file2)
    }) %...>% (function(filtered_files) {
      # Update reactive values with the newly filtered files
      filtered_files(list(file1 = filtered_files$filtered_file1, file2 = filtered_files$filtered_file2))
      
      # Notify the user that filtering has completed
      showNotification("Filtering reads below quality threshold completed.", duration = 2)
    }) %...!% (function(err) {
      # Notify the user if an error occurs during filtering
      showModal(modalDialog(
        title = "Error",
        paste("Error during filtering process:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
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
    phred_offset <- get_phred_offset(input$quality_threshold)
    
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
    req(input$genome_fasta_indexing, input$memory_size, input$index_base_name)
    
    genome_fasta <- input$genome_fasta_indexing$datapath
    memory_size <- input$memory_size
    index_base_name <- input$index_base_name  # Use the basename provided by the user
    working_directory <- working_dir()
    
    showNotification("Starting genome index creation. This may take between 1.5 to 2 hours.", duration = NULL)
    
    future({
      setwd(working_directory)
      
      # Build the index with the basename chosen by the user
      buildindex(
        basename = index_base_name,
        reference = genome_fasta,
        gappedIndex = FALSE,
        indexSplit = TRUE,
        memory = memory_size,
        TH_subread = 100,
        colorspace = FALSE
      )
      
      # Get the list of index files created
      index_files <- list.files(pattern = paste0(index_base_name, ".*\\.(reads|files|txt|b.array|b.tab)"), 
                                path = working_directory, full.names = TRUE)
      
      list(index_files = index_files, index_base_name = index_base_name)
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
  
  observeEvent(input$filter_reads, {
    # Ensure either demultiplexed or original FASTQs are available
    req(
      (isTruthy(demultiplexed_files()$file1) && isTruthy(demultiplexed_files()$file2)) ||
        (isTruthy(fastq_files()$file1) && isTruthy(fastq_files()$file2))
    )
    
    # Notify the user that filtering has started
    showNotification("Filtering reads below quality threshold started.", duration = 2)
    
    # Retrieve the phred_offset and quality_threshold inside the reactive context
    phred_offset <- get_phred_offset(input$machine_type)
    quality_threshold <- input$quality_threshold
    
    # Define the output directory
    output_dir <- working_dir()
    
    # Determine which files to use: demultiplexed or original
    file_to_use1 <- if (isTruthy(demultiplexed_files()$file1)) {
      demultiplexed_files()$file1$datapath
    } else {
      fastq_files()$file1$datapath
    }
    
    file_to_use2 <- if (isTruthy(demultiplexed_files()$file2)) {
      demultiplexed_files()$file2$datapath
    } else {
      fastq_files()$file2$datapath
    }
    
    future({
      # Inside the future, we now have access to the values of phred_offset and quality_threshold
      
      # Filter the reads for each FASTQ file with a timestamped name
      filtered_file1 <- filter_fastq_by_quality(file_to_use1, output_dir, phred_offset, quality_threshold)
      filtered_file2 <- filter_fastq_by_quality(file_to_use2, output_dir, phred_offset, quality_threshold)
      
      list(filtered_file1 = filtered_file1, filtered_file2 = filtered_file2)
    }, seed = TRUE) %...>% (function(filtered_files) {
      # Update reactive values with the newly filtered files
      filtered_files(list(file1 = filtered_files$filtered_file1, file2 = filtered_files$filtered_file2))
      
      # Notify the user that filtering has completed
      showNotification("Filtering reads below quality threshold completed.", duration = 2)
    }) %...!% (function(err) {
      # Notify the user if an error occurs during filtering
      showModal(modalDialog(
        title = "Error",
        paste("Error during filtering process:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  

  
  observeEvent(input$demultiplex, {
    req(input$barcodes_file, input$fastq_file1)
    
    # Notify the user that demultiplexing has started
    showNotification("Demultiplexing process started.", duration = 2)
    
    # Retrieve the phred_offset and pass it to future
    phred_offset <- get_phred_offset(input$machine_type)
    
    barcodes_file <- input$barcodes_file$datapath
    fastq_file1 <- input$fastq_file1$datapath
    output_dir <- working_dir()
    
    future({
      # Inside the future, we now have access to the value of phred_offset
      barcodes <- read.csv(barcodes_file, header = FALSE, stringsAsFactors = FALSE, skip = 1)
      barcodes <- as.character(barcodes[[1]])
      barcodes <- gsub("\\s", "", barcodes)
      
      barcode_lengths <- unique(nchar(barcodes))
      if (length(barcode_lengths) > 1) {
        stop("All barcodes must be of the same length.")
      }
      
      # Use the phred_offset value directly
      custom_demultiplex(fastq_file1, barcodes, output_dir, phred_offset)
      
      TRUE
    }, seed = TRUE) %...>% (function(result) {
      showNotification("Demultiplexing completed.", duration = 2)
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
    
    # Notify the user that GFF3 generation has started
    showNotification("Generating GFF3 file started.", duration = 2)
    
    gtf_file <- input$gtf_file$datapath
    working_directory <- working_dir()
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    future({
      # Generate the GFF3 file
      gff_file <- file.path(working_directory, paste0("annotations_", timestamp, ".gff3"))
      gtf <- import(gtf_file, format = "gtf")
      export(gtf, gff_file, format = "gff3")
      gff_file
    }, seed = TRUE) %...>% (function(gff_file) {
      # Notify the user that GFF3 generation was successful
      showNotification("GFF3 file generated successfully.", duration = 2)
      
      # Display the path to the generated GFF3 file
      output$gff3_file_path <- renderText({
        paste("Generated GFF3 File Path: ", gff_file)
      })
    }) %...!% (function(err) {
      # Notify the user if an error occurs during GFF3 generation
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
