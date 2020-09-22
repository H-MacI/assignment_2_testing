#Bash script designed to organize the files contained in the sample folder given in the assignment#
mkdir src
mkdir src/pahfit_programs
mkdir src/python_programs

git mv *ipynb src/python_programs
git mv *.pro src/pahfit_programs
#################################
mkdir data
mkdir data/fits_files
mkdir data/dat_files
mkdir data/excel_files
mkdir data/ds9_files
mkdir data/tbl_files
mkdir data/incorrect
mkdir data/txt_files

git mv *.dat data/dat_files
git mv *.fits data/fits_files
git mv *.xls data/excel_files
git mv *.tbl data/tbl_files
git mv *.reg data/ds9_files
git mv *Wrong/*.tbl data/incorrect
git mv *Wrong/*.xdr data/incorrect
git mv nucUNC data/txt_files
git mv nucFLUX data/txt_files
git mv *phot.txt data/txt_files
git mv *phot data/txt_files
#################################
mkdir docs
mkdir docs/session_logs
mkdir docs/analysis
mkdir docs/file_types

git mv *.log docs/session_logs
git mv *journal.txt docs/session_logs
git mv analysis*.txt docs/analysis
git mv IRS*2015.txt docs/analysis
##################################
mkdir results
mkdir results/images
mkdir results/images/eps_files
mkdir results/images/png_files
mkdir results/images/ps_files
mkdir results/images/incorrect
mkdir results/fitting_parameters
mkdir results/fitting_parameters/correct
mkdir results/fitting_parameters/incorrect
mkdir results/structures
mkdir results/tbl_files


git mv *report.txt results/fitting_parameters/incorrect
git mv *correct.txt results/fitting_parameters/correct
git mv *pahfit.png results/images
git mv *.xdr results/structures
git mv *.png results/images/png_files
git mv *.eps results/images/eps_files
git mv *.ps results/images/ps_files
git mv *Wrong/*.png results/images/incorrect
#####################################
if [ -e README -a -e README2 ]
then
    cat README README2 > File_Contents
    git rm README
    git rm README2
fi
git add File_Contents
git mv *File_Contents docs/file_types

