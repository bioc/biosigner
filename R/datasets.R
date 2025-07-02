#' Analysis of plasma from diabetic patients by LC-HRMS
#'
#' Plasma samples from 69 diabetic patients were analyzed by reversed-phase
#' liquid chromatography coupled to high-resolution mass spectrometry (Orbitrap
#' Exactive) in the negative ionization mode. The raw data were pre-processed
#' with XCMS and CAMERA (5,501 features), corrected for signal drift, log10
#' transformed, and annotated with an in-house spectral database. The patient's
#' age, body mass index, and diabetic type are recorded. These three clinical
#' covariates are strongly associated, most of the type II patients being older
#' and with a higher bmi than the type I individuals.
#'
#' @name diaplasma
#' @docType data
#' @format A list with the following elements:
#' \describe{
#' \item{dataMatrix}{a 69 samples x
#' 5,501 features matrix of numeric type corresponding to the intensity
#' profiles (values have been log10-transformed),}
#' \item{sampleMetadata}{a 69 x 3 data frame, with the patients' diabetic type
#' ('type', factor), age ('age', numeric), and body mass index ('bmi', numeric),}
#' \item{variableMetadata}{a 5,501 x 8 data frame, with the median m/z ('mzmed',
#' numeric) and the median retention time in seconds ('rtmed', numeric) from
#' XCMS, the 'isotopes' (character), 'adduct' (character) and 'pcgroups'
#' (numeric) annotations from CAMERA, and the names of the m/z and RT matching
#' compounds from an in-house database of pure spectra from commercial
#' metabolites ('spiDb', character).}
#' }
#' @return List containing the 'dataMatrix' matrix (numeric) of data (samples
#' as rows, variables as columns), the 'sampleMetadata' data frame of sample
#' metadata, and the variableMetadata data frame of variable metadata. Row
#' names of 'dataMatrix' and 'sampleMetadata' are identical. Column names of
#' 'dataMatrix' are identical to row names of 'variableMetadata'. For details
#' see the 'Format' section above.
#' @references Rinaudo P., Boudah S., Junot C. and Thevenot E.A. (2016).
#' biosigner: a new method for the discovery of significant molecular signatures
#' from omics data. Frontiers in Molecular Biosciences 3.
#' doi:10.3389/fmolb.2016.00026
#' @source 'diaplasma' dataset.
#' @keywords datasets
NULL

#' Spike-in metabolomics data for apple extracts (from the BioMark package)
#'
#' Data from a spike-in experiment for apple extracts. Twenty apple extracts are
#' divided in two groups, one control, and one spike-in group. The control group
#' is measured without any spiking - the spike-in group is spiked with nine
#' chemical compounds in three different combinations of concentrations.
#' The data provide the experimental data of the forty apple extracts in the
#' list 'SpikePos' for positive ionization, and in the separate data.frame 'pos.markers'
#' that contains information about the features of the standards, i.e., the spike-in compounds.
#' The data use CAMERA grouping to automatically determine which features are
#' corresponding to which spike-in compounds. Raw data in CDF format are available
#' from the MetaboLights repository (MTBLS59).
#'
#' @name SpikePos
#' @docType data
#' @format 'SpikePos' is a list with the following three elements:
#' \describe{
#' \item{data}{a 40 samples x 1,632 features matrix of numeric type, describing for each of the forty injections the intensity of the features (columns). Column names consist of a combination of retention time (in seconds) and m/z values, and are sorted on retention time,}
#' \item{classes}{a factor containing the class labels for the forty injections (control, or group1, 2 or 3),}
#' \item{annotation}{a 1,632 features x 11 metadata data.frame, containing for each of the features XCMS and CAMERA information, such as mz, rt, number of times a feature is identified in the control or spike-in samples, possible isotope or adduct annotation, and whether or not the feature is identified in the standards (the spike-in data).}
#' }
#' In addition, 'pos.markers' is a data frame that contains the information of the standards, i.e. the compounds that are spiked in. These data.frames describe in their rows single features identified with XCMS and CAMERA, using the same settings as the experimental apple data, and have the following columns:
#' \describe{
#' \item{comp}{The (short) name of the spiked-in compound giving rise to this particular feature,}
#' \item{mz, rt, isotope, adduct}{Feature information, similar to the information in the 'annotation' fields in 'SpikePos',}
#' \item{feature.nr}{The number of the corresponding feature in 'SpikePos',}
#' \item{group1, group2, group3}{Approximate spiking levels for the three groups. A value of 1.0 corresponds to an increase that is roughly equal to the naturally occuring concentration in apple. Exceptions are trans-resveratrol and cyanidin-3-galactoside, both not naturally occuring. These two compounds have been spiked in at one constant level which gives features of comparable size.}
#' }
#' @return 'SpikePos' list and 'pos.markers' data frame. For details see the 'Format' section above.
#' @examples
#' data(SpikePos)
#' plot(SpikePos$annotation[,c('rt', 'mz')], xlab = 'Time (s)', ylab = 'm/z', main = 'Positive ionization mode')
#' points(pos.markers[!is.na(pos.markers$feature.nr), c('rt', 'mz')], pch = 19, col = 2)
#' @author Pietro Franceschi
#' @references Franceschi,P. et al. (2012) A benchmark spike-in data set for biomarker identification in metabolomics. Journal of Chemometrics, 26, 16–24. DOI:10.1002/cem.1420.
#' @source 'SpikePos' dataset.
#' @keywords datasets
NULL

