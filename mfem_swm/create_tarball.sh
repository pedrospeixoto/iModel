#! /bin/bash
# NEEDS TO BE UPDATED!! PXT
PACKAGE_DIRS_AND_FILES="
doc
nml
gmt
matlab
grd/*gmt
"

PACKAGE_DIRS_AND_FILES="$PACKAGE_DIRS_AND_FILES 
dodecahedron.xref
adjacency_matrices.f90
buildop_fem.f90
femswe.f90
gengrid_cube.f90
gengrid_hex.f90
sfc_optimize.f90
findevals2.m  jtaxes.m     jtplotgrid.m  plotdiag.m
join.m        jtcontour.m  jtrotplot.m
Makefile
create_all.sh
grid.f90
compare_femswe_hex.sh
compare_femswe_cube.sh
create_jobs.py
computeHilbert2d.c
computeHilbert3d.c
"
 

PACKAGE_NAME="pdfemcode3_`date +%Y_%m_%d`"
PACKAGE_DIR="$PACKAGE_NAME"
PACKAGE_TARBALL="$PACKAGE_NAME.tar.bz2"

echo "Creating package $PACKAGE_NAME"
rm -f -r "$PACKAGE_DIR"
mkdir "$PACKAGE_DIR"

echo " + copying files"
for file in $PACKAGE_DIRS_AND_FILES; do
	cp -r "../$file" "$@" "$PACKAGE_DIR"
done

echo " + removing svn information"
# remove svn from package directory
cd "$PACKAGE_DIR" && { find ./ -name ".svn" | xargs rm -Rf; } && cd ..

echo " + creating tarball $PACKAGE_TARBALL"
rm -f "$PACKAGE_TARBALL"
tar cjf "$PACKAGE_TARBALL" "$PACKAGE_DIR"

echo " + cleaning up"
rm -r "$PACKAGE_DIR"

echo "Done!" 
