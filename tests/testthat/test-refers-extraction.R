context("reference extraction")
test_that("malformed input error",{
    td <- tempdir()
    td <- gsub("[\\]","/",td)
    refs <- dir(system.file(package="Rhisat", "extdata", "bt2","refs"),full=TRUE)
    idx <- file.path(td, "lambda_virus")
    reads_1 <- system.file(package="Rhisat", "extdata", "bt2", "reads","reads_1.fastq")
    reads_2 <- system.file(package="Rhisat", "extdata", "bt2", "reads","reads_2.fastq")

    hisat_build(references=refs, bt2Index=idx,"--quiet",overwrite=TRUE)

    expect_error(
        hisat_inspect(bt2Index=TRUE, "--summary")
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
        hisat_inspect(
            bt2Index = idx,"--summary", exe=FALSE),
        paste(
            file.path(system.file(package="Rhisat"), "hisat-inspect"),
            "--summary", idx
        )
    )
}
)