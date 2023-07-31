#!/bin/bash
cp Generic_Grains2D.dtd Grains2D.dtd
cp Generic_Grains3D.dtd Grains3D.dtd
cp Generic_Grains2D_InFluid.dtd Grains2D_InFluid.dtd
cp Generic_Grains3D_InFluid.dtd Grains3D_InFluid.dtd

contact=${GRAINS_HOME}"/Main/dtd/Convex.dtd"
comm="s|___WHERESISCONVEX__|$contact|g"

sed -e $comm Grains2D.dtd > __tttt__
mv __tttt__ Grains2D.dtd
sed -e $comm Grains3D.dtd > __tttt__
mv __tttt__ Grains3D.dtd
sed -e $comm Grains2D_InFluid.dtd > __tttt__
mv __tttt__ Grains2D_InFluid.dtd
sed -e $comm Grains3D_InFluid.dtd > __tttt__
mv __tttt__ Grains3D_InFluid.dtd

contact=${GRAINS_HOME}"/Main/dtd/Contact.dtd"
comm="s|___WHERESISCONTACT__|$contact|g"

sed -e $comm Grains2D.dtd > __tttt__
mv __tttt__ Grains2D.dtd
sed -e $comm Grains3D.dtd > __tttt__
mv __tttt__ Grains3D.dtd
sed -e $comm Grains2D_InFluid.dtd > __tttt__
mv __tttt__ Grains2D_InFluid.dtd
sed -e $comm Grains3D_InFluid.dtd > __tttt__
mv __tttt__ Grains3D_InFluid.dtd
