context("read alignment")
test_that("malformed input error",{
    td <- tempdir()
    td <- gsub("[\\]","/",td)
    refs <- dir(system.file(package="Rhisat", "extdata", "bt2","refs"),full=TRUE)
    idx <- file.path(td, "lambda_virus")
    reads_1 <- system.file(package="Rhisat", "extdata", "bt2", "reads","reads_1.fastq")
    reads_2 <- system.file(package="Rhisat", "extdata", "bt2", "reads","reads_2.fastq")

    ## Avoid warning of overwrite
    options (warn = -1)
    hisat_build(references=refs, bt2Index=idx,"--quiet",overwrite=TRUE)

    expect_error(
        hisat(bt2Index = TRUE,
        samOutput = file.path(td, "result.sam"),
        seq1=reads_1,seq2=reads_2,overwrite=TRUE,"--threads 3")
    )

    expect_error(
        hisat(bt2Index = idx,
        samOutput = NULL,
        seq1=reads_1,seq2=reads_2,overwrite=TRUE,"--threads 3")
    )

    expect_error(
        hisat(bt2Index = idx,
        samOutput = file.path(td, "result.sam"),
        seq1=NULL,seq2=reads_2,overwrite=TRUE,"--threads 3")
    )

    expect_error(
        hisat(bt2Index = idx,
        samOutput = file.path(td, "result.sam"),
        seq1=reads_1,seq2=reads_2,"-U")
    )

}
)

test_that("formatted input work",{
    td <- tempdir()
    td <- gsub("[\\]","/",td)
    refs <- dir(system.file(package="Rhisat", "extdata", "bt2","refs"),full=TRUE)
    idx <- file.path(td, "lambda_virus")
    reads_1 <- system.file(package="Rhisat", "extdata", "bt2", "reads","reads_1.fastq")
    reads_2 <- system.file(package="Rhisat", "extdata", "bt2", "reads","reads_2.fastq")
    hisat_build(references=refs, bt2Index=idx,"--quiet",overwrite=TRUE)

    expect_equal(
        hisat(
            bt2Index = idx,
            samOutput = file.path(td, "result.sam"),
            seq1=reads_1,seq2=reads_2,overwrite=TRUE,
            exe=FALSE),
        paste(
            file.path(system.file(package="Rhisat"), "hisat"),
            "-x", idx,
            "-1", reads_1,
            "-2", reads_2,
            "-S", file.path(td, "result.sam")
        )
    )

    expect_equal(
        hisat(
            bt2Index = idx,
            samOutput = file.path(td, "result.sam"),
            seq1=reads_1,seq2=NULL,overwrite=TRUE,
            exe=FALSE),
        paste(
            file.path(system.file(package="Rhisat"), "hisat"),
            "-x", idx,
            "-U", reads_1,
            "-S", file.path(td, "result.sam")
        )
    )

    cmdout <- hisat(bt2Index = idx,samOutput = file.path(td, "result.sam"),seq1=reads_1,overwrite=TRUE)
    expect_is(cmdout, "character")
    res <- head(readLines(file.path(td, "result.sam")))
    expect_length(res, 6)

}
)