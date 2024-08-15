# Load required libraries
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

# Custom demultiplexing function
custom_demultiplex <- function(fastq_file, barcodes, output_dir) {
  sequences <- sread(readFastq(fastq_file))
  qualities <- extract_quality_scores(fastq_file)
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
  
  for (barcode in names(demux_results)) {
    output_file <- file.path(output_dir, paste0("demultiplexed_", barcode, ".fastq"))
    seqs <- demux_results[[barcode]]$sequences
    quals <- demux_results[[barcode]]$qualities
    
    quality_scores <- BStringSet(sapply(quals, function(q) {
      paste(convert_numeric_to_ascii(q), collapse = "")
    }))
    writeFastq(ShortReadQ(sread = DNAStringSet(seqs), quality = quality_scores), output_file, mode = "w", compress = FALSE, append = FALSE)
  }
  
  return(output_dir)
}

# Function to read a FASTQ file and extract quality scores
extract_quality_scores <- function(fastq_file) {
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
      as.integer(charToRaw(char)) - 33
    })
  }
  
  quality_scores <- lapply(quality_scores, convert_to_phred)
  return(quality_scores)
}

# Function to convert numeric Phred scores to ASCII characters
convert_numeric_to_ascii <- function(numeric_qualities) {
  sapply(numeric_qualities, function(score) {
    rawToChar(as.raw(score + 33))
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


ui <- fluidPage(
  theme = shinytheme("cyborg"),
  titlePanel("JCap Paired_End Fastq Processing App"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("directory", "Set Working Directory", "Please select a directory"),
      fileInput("fastq_file1", "Upload FASTQ File 1", accept = c(".fastq", ".fq")),
      fileInput("fastq_file2", "Upload FASTQ File 2", accept = c(".fastq", ".fq")),
      textInput("sample_name", "Name of Sample", value = "Sample1"),
      fileInput("barcodes_file", "Upload Barcodes File (Optional)", accept = c(".csv", ".tsv", ".txt")),
      actionButton("demultiplex", "Run Demultiplexing"),
      fileInput("demultiplexed_file1", "Upload Demultiplexed File 1", accept = c(".fastq", ".fq")),
      fileInput("demultiplexed_file2", "Upload Demultiplexed File 2", accept = c(".fastq", ".fq")),
      selectInput("machine_type", "Select Sequencing Machine Type", 
                  choices = c("Illumina (Phred+33)" = "phred33", 
                              "Sanger (Phred+33)" = "phred33", 
                              "Solexa/Illumina (Phred+64)" = "phred64")),
      actionButton("generate_summary", "Generate FASTQ Summary"),
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


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 2 * 1024^50000)
  
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
  
  observeEvent(input$generate_summary, {
    files_to_use <- if (!is.null(filtered_files()$file1) && !is.null(filtered_files()$file2)) {
      c(filtered_files()$file1$datapath, filtered_files()$file2$datapath)
    } else if (!is.null(demultiplexed_files()$file1) && !is.null(demultiplexed_files()$file2)) {
      c(demultiplexed_files()$file1$datapath, demultiplexed_files()$file2$datapath)
    } else {
      c(fastq_files()$file1$datapath, fastq_files()$file2$datapath)
    }
    req(files_to_use)
    showNotification("Generating FASTQ summary", duration = 2)
    
    filename1 <- basename(files_to_use[1])
    filename2 <- if (length(files_to_use) > 1) basename(files_to_use[2]) else NA
    quality_threshold <- input$quality_threshold
    phred_offset <- get_phred_offset(input$machine_type)
    
    future({
      quality_scores1 <- extract_quality_scores(files_to_use[1])
      if (!is.na(filename2)) {
        quality_scores2 <- extract_quality_scores(files_to_use[2])
      } else {
        quality_scores2 <- NULL
      }
      
      res1 <- process_chunk(quality_scores1, quality_threshold)
      if (!is.null(quality_scores2)) {
        res2 <- process_chunk(quality_scores2, quality_threshold)
      } else {
        res2 <- NULL
      }
      
      avg_quality1 <- res1$avg_quality
      percent_above_threshold1 <- res1$percent_above_threshold
      num_reads1 <- res1$num_reads
      
      if (!is.null(res2)) {
        avg_quality2 <- res2$avg_quality
        percent_above_threshold2 <- res2$percent_above_threshold
        num_reads2 <- res2$num_reads
      } else {
        avg_quality2 <- NA
        percent_above_threshold2 <- NA
        num_reads2 <- NA
      }
      
      summary_df <- data.frame(
        Filename = c(filename1, if (!is.na(filename2)) filename2 else NULL),
        NumReads = c(num_reads1, if (!is.null(num_reads2)) num_reads2 else NULL),
        PercentAboveThreshold = c(percent_above_threshold1, if (!is.null(percent_above_threshold2)) percent_above_threshold2 else NULL),
        AvgQuality = c(avg_quality1, if (!is.null(avg_quality2)) avg_quality2 else NULL)
      )
      summary_df
    }, seed = TRUE) %...>% (function(summary_df) {
      output$fastq_summary <- renderDT({
        datatable(summary_df)
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
  
  observeEvent(input$filter_reads, {
    # Capture the current value of working_dir() in a non-reactive variable
    current_working_dir <- working_dir()
    
    # Determine which files to use for filtering
    files_to_use <- if (!is.null(filtered_files()$file1) && !is.null(filtered_files()$file2)) {
      c(filtered_files()$file1$datapath, filtered_files()$file2$datapath)
    } else if (!is.null(demultiplexed_files()$file1) && !is.null(demultiplexed_files()$file2)) {
      c(demultiplexed_files()$file1$datapath, demultiplexed_files()$file2$datapath)
    } else {
      c(fastq_files()$file1$datapath, fastq_files()$file2$datapath)
    }
    req(files_to_use)
    req(input$machine_type)
    
    showNotification("Filtering reads", duration = 2)
    
    quality_threshold <- input$quality_threshold
    
    future({
      # Get the current date and time for unique filenames
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      # Iterate over the list of files, filter each one, and write the output FASTQ file
      filtered_files_list <- lapply(seq_along(files_to_use), function(i) {
        output_filename <- paste0(timestamp, "_filtered_", i, ".fastq")
        filter_and_write_fastq(files_to_use[[i]], quality_threshold, current_working_dir, output_filename)
      })
      
      filtered_files_list
    }, seed = TRUE) %...>% (function(result) {
      showNotification("Reads filtered successfully.")
      
      # Get read counts for each file
      read_counts <- sapply(result, get_read_count)
      
      # Update the filtered_files reactive list with the new filtered files
      filtered_files(list(file1 = result[[1]], file2 = ifelse(length(result) > 1, result[[2]], NULL)))
      
      # Notify the user about the number of reads in each file and their paths
      showNotification(paste("Filtered FASTQ 1 has", read_counts[1], "reads. Path:", result[[1]]), duration = 5)
      if (length(result) > 1 && !is.null(result[[2]])) {
        showNotification(paste("Filtered FASTQ 2 has", read_counts[2], "reads. Path:", result[[2]]), duration = 5)
      }
    }) %...!% (function(err) {
      showModal(modalDialog(
        title = "Error",
        paste("Error filtering reads:", err$message),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  # Function to filter reads, write them to a new FASTQ file, and return the file path
  filter_and_write_fastq <- function(fastq_file, quality_threshold, output_dir, output_filename) {
    phred_offset <- 33
    reads <- readFastq(fastq_file)
    quality_scores <- extract_quality_scores(fastq_file)
    
    # Filter reads based on the quality threshold
    pass_filter <- sapply(quality_scores, function(q) mean(q) >= quality_threshold)
    filtered_reads <- reads[pass_filter]
    
    # Construct the full path for the filtered file in the working directory
    filtered_file <- file.path(output_dir, output_filename)
    
    # Write the filtered reads to the new file
    writeFastq(filtered_reads, filtered_file, mode = "w", compress = FALSE, append = FALSE)
    
    # Return the path of the filtered file
    return(filtered_file)
  }
  
  # Function to get the number of reads in a FASTQ file
  get_read_count <- function(fastq_file) {
    reads <- readFastq(fastq_file)
    return(length(reads))
  }
  
  observeEvent(input$create_index, {
    req(input$genome_fasta_indexing, input$memory_size)
    
    genome_fasta <- input$genome_fasta_indexing$datapath
    memory_size <- input$memory_size
    working_directory <- working_dir()
    
    showNotification("Starting genome index creation. This may take between 1.5 to 2 hours.", duration = NULL)
    
    future({
      setwd(working_directory)
      
      # Use timestamp for a unique base name
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
    # Capture the current value of reactive inputs inside the reactive context
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
          isPairedEnd = TRUE  # Always using paired-end reads
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
      message("Starting demultiplexing")
      
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
