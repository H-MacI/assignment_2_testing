mkdir src
mkdir src/pahfit_programs
mkdir data
mkdir data/fits_files
mkdir data/dat_files
mkdir data/excel_files
mkdir data/ds9_files
mkdir data/tbl_files
mkdir data/incorrect

mkdir docs

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
git mv *.pro src/pahfit_programs

git mv *.dat data/dat_files
git mv *.fits data/fits_files
git mv *.xls data/excel_files
git mv *.tbl results/tbl_files
git mv *.reg data/ds9_files
git mv *Wrong/*.tbl data/incorrect
git mv *Wrong/*.xdr data/incorrect

