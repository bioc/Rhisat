#' @name hisat
#' @title Interface to hisat of hisat-0.1.6-beta-source.zip
#' @description This function can be use to call wrapped \code{hisat}
#' binary.
#' @param bt2Index \code{Character} scalar. hisat index files
#' prefix: 'dir/basename'
#' (minus trailing '.*.bt2' of 'dir/basename.*.bt2').
#' @param samOutput \code{Character} scalar. A path to a SAM file
#' used for the alignment output.
#' @param seq1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates
#' paired with file paths in seq2.
#' @param seq2 \code{Character} vector. It contains file paths with
#' #2 mates paired with file paths in seq1.
#' @param ... Additional arguments to be passed on to the binaries.
#' See below for details.
#' @param overwrite \code{Logical}. Force overwriting of existing
#' files if setting \code{TRUE}.
#' @details All additional arguments in ... are interpreted as
#' additional parameters to be passed on to
#' hisat. All of them should be \code{Character} or
#' \code{Numeric} scalar. You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8 --no-mixed")
#' with white space splited just like command line,
#' or put them in different \code{Character}
#' (e.g. "--threads","8","--no-mixed"). Note that some
#' arguments("-x","-U","-1","-2","-S") to the
#' hisat are invalid if they are already handled as explicit
#' function arguments. See the output of
#' \code{hisat_usage()} for details about available parameters.
#' @author QiXiu Du
#' @return An invisible \code{Integer} of call
#' status. The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012).
#' Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export hisat
#' @examples
#' td <- tempdir()
#' ## Building a hisat index
#' refs <- dir(system.file(package="Rhisat", "extdata", "bt2","refs"),
#' full=TRUE)
#' hisat_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads 4 --quiet",overwrite=TRUE)
#' ## Alignments
#' reads_1 <- system.file(package="Rhisat", "extdata", "bt2", "reads",
#' "reads_1.fastq")
#' reads_2 <- system.file(package="Rhisat", "extdata", "bt2", "reads",
#' "reads_2.fastq")
#' if(file.exists(file.path(td, "lambda_virus.1.bt2"))){
#'     cmdout<-hisat(bt2Index = file.path(td, "lambda_virus"),
#'        samOutput = file.path(td, "result.sam"),
#'        seq1=reads_1,seq2=reads_2,overwrite=TRUE,"--threads 3");cmdout
#'     head(readLines(file.path(td, "result.sam")))
#' }
#'

hisat <- function(bt2Index,samOutput,seq1,...,seq2=NULL,overwrite=FALSE){
    if(R.Version()$arch=="i386"){
        return("hisat is not available for 32bit, please use 64bit R instead")
    }
    bt2Index <-trimws(as.character(bt2Index))
    samOutput<-trimws(as.character(samOutput))

    seq1<-trimws(as.character(seq1))


    if(!is.null(seq2)){
        seq2<-trimws(as.character(seq2))
        if(length(seq1)!=length(seq2)){
            stop(paste0("The lengths of arguments ",
                        "`seq1` and `seq2` should be the same length"))
        }
    }
    paramArray<-checkAddArgus("-x|-U|-1|-2|-S",...)


    checkFileExist(seq1,"seq1")
    checkFileExist(seq2,"seq2")
    checkPathExist(bt2Index,"bt2Index")
    checkFileExist(paste0(bt2Index,".1.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".2.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".3.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".4.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".rev.1.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".rev.2.bt2"),"bt2Index")
    checkFileCreatable(samOutput,"samOutput",overwrite)

    argvs = c("-x",bt2Index)
    seq1<-paste0(seq1,collapse = ",")
    if(is.null(seq2)){
        argvs <- c(argvs,"-U",seq1)
    }else{
        seq2<-paste0(seq2,collapse = ",")
        argvs <- c(argvs,"-1",seq1,"-2",seq2)
    }

    argvs <- c(paramArray,argvs,"-S",samOutput)

    invisible(.callbinary("hisat-align-s",paste(argvs,collapse = " ")))

}


#' @name hisat_build
#' @title Interface to hisat-build of hisat-0.1.6-beta
#' @description This function can be use to call wrapped \code{hisat-build}
#' binary
#' @param references \code{Character} vector. The path to the files containing
#' the references for which to
#' build a bowtie index.
#' @param bt2Index \code{Character} scalar. Write hisat index data to files
#' with this prefix: 'dir/basename'.
#' If the files with path like 'dir/basename.*.bt2' already exists,
#' the function function will cast an error,
#' unless argument overwrite is \code{TRUE}.
#' @param ... Additional arguments to be passed on to the binaries.
#' See below for details.
#' @param overwrite \code{Logical}. Force overwriting of existing files
#' if setting \code{TRUE}.
#' @details All additional arguments in ... are interpreted as additional
#' parameters to be passed on to
#' hisat_build. All of them should be \code{Character} or
#' \code{Numeric} scalar. You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8 --quiet") with white
#' space splited just like command line,
#' or put them in different \code{Character}(e.g. "--threads","8","--quiet").
#' See the output of
#' \code{hisat_build_usage()} for details about available parameters.
#' @author Qixiu Du
#' @return An invisible \code{Integer} of call status.
#' The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read
#' alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export hisat_build
#' @examples
#' td <- tempdir()
#' ## Building a hisat index
#' refs <- dir(system.file(package="Rhisat", "extdata", "bt2","refs"),
#' full=TRUE)
#' cmdout<-hisat_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads 4 --quiet",overwrite=TRUE);cmdout
#' ## Use additional arguments in another way
#' cmdout<-hisat_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads",4,"--quiet",overwrite=TRUE);cmdout
#' ## The function will print the output
#' ## during the process without "--quiet" argument.
#' cmdout<-hisat_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' overwrite=TRUE);cmdout

hisat_build <- function(references,bt2Index,...,overwrite=FALSE){
    if(R.Version()$arch=="i386"){
        return("hisat is not available for 32bit, please use 64bit R instead")
    }
    references<- trimws(as.character(references))
    bt2Index <- trimws(as.character(bt2Index))

    paramArray<-checkAddArgus("noinvalid",...)

    checkFileExist(references,"references")
    checkPathExist(bt2Index,"bt2Index")
    checkFileCreatable(paste0(bt2Index,".1.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".2.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".3.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".4.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".rev.1.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".rev.2.bt2"),"bt2Index",overwrite)

    references<-paste0(references,collapse = ",")
    argvs <- c(paramArray,references,bt2Index)


    invisible(.callbinary("hisat-build-s",paste(argvs,collapse = " ")))

}

#' @name hisat_version
#' @title Print version information of hisat-0.1.6-beta
#' @description Print version information of hisat-0.1.6-beta
#' @author Qixiu Du
#' @return An invisible \code{Integer} of call status.
#' The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read
#' alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export hisat_version
#' @examples
#' cmdout<-hisat_version();cmdout
hisat_version <- function(){
    if(R.Version()$arch=="i386"){
        return("hisat is not available for 32bit, please use 64bit R instead")
    }
    .callbinary("hisat-align-s","--version")
}

#' @name hisat_usage
#' @title Print available arguments for hisat
#' @description Note that some arguments to the
#' hisat are invalid if they are
#' already handled as explicit function arguments.
#' @author Qixiu Du
#' @return hisat available arguments and their usage.
#' @references Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read
#' alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export hisat_usage
#' @examples
#' hisat_usage()
hisat_usage <- function(){
    if(R.Version()$arch=="i386"){
        return("hisat is not available for 32bit, please use 64bit R instead")
    }
    .callbinary("hisat-align-s","-h")
}

#' @name hisat_build_usage
#' @title Print available arguments for hisat_build_usage
#' @description Note that some arguments to the
#' hisat_build_usage are invalid if they are
#' already handled as explicit function arguments.
#' @author Qixiu Du
#' @return hisat_build available arguments and their usage.
#' @references Langmead B, Salzberg S.
#' Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
#' @export hisat_build_usage
#' @examples
#' hisat_build_usage()
hisat_build_usage <- function() {
    if(R.Version()$arch=="i386"){
        return("hisat is not available for 32bit, please use 64bit R instead")
    }
    .callbinary("hisat-build-s","-h")
}


#' @name hisat_inspect
#' @title Interface to hisat-inspect of hisat-0.1.6-beta
#' @description This function can be use to call wrapped \code{hisat-inspect}
#' binary.
#' @param bt2_base \code{Character} scalar.name of any of the index files 
#' but with the .X.bt2 or .rev.X.bt2 suffix omitted.
#' @param ... Additional arguments to be passed on to the binaries.
#' See below for details.
#' @author QiXiu Du
#' @return a FASTA file containing the sequences of the original references 
#' (with all non-A/C/G/T characters converted to Ns)
#' @references Langmead, B., & Salzberg, S. L. (2012).
#' Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export hisat_inspect
#' @examples
#' td <- tempdir()
#' ## Building a hisat index
#' refs <- dir(system.file(package="Rhisat", "extdata", "bt2","refs"),
#' full=TRUE)
#' hisat_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads 4 --quiet",overwrite=TRUE)
#' if(file.exists(file.path(td, "lambda_virus.1.bt2"))){
#'     cmdout<-hisat_inspect(bt2Index = file.path(td, "lambda_virus"),
#'        "--summary");cmdout 
#' }
#'

hisat_inspect <- function(bt2Index,...){
    if(R.Version()$arch=="i386"){
        return("hisat is not available for 32bit, please use 64bit R instead")
    }
    bt2Index <-trimws(as.character(bt2Index))
    paramArray<-checkAddArgus("noinvalid",...)
    checkFileExist(paste0(bt2Index,".1.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".2.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".3.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".4.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".rev.1.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".rev.2.bt2"),"bt2Index")

    argvs <- c(paramArray,bt2Index)

    invisible(.callbinary("hisat-inspect-s",paste(argvs,collapse = " ")))

}