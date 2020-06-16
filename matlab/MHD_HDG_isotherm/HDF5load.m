%****************************************************************************
%  Copyright Euratom-CEA
%  Authors : Virginie Grandgirard (virginie.grandgirard@cea.fr)
%            Chantal Passeron (chantal.passeron@cea.fr)
%            Guillaume Latu (latu.labri@cea.fr)
%            Xavier Garbet (xavier.garbet@cea.fr)
%            Philippe Ghendrih (philippe.ghendrih@cea.fr)
%            Yanick Sarazin (yanick.sarazin@cea.fr)
%
%  This code GYSELA (for GYrokinetic SEmi-LAgrangian) is a 5D
%  gyrokinetic global full-f code for simulating the plasma turbulence
%  in a tokamak.
%
%  This software is governed by the CeCILL-B license under French law and
%  abiding by the rules of distribution of free software.  You can  use,
%  modify and redistribute the software under the terms of the CeCILL-B
%  license as circulated by CEA, CNRS and INRIA at the following URL
%  "http://www.cecill.info".
%
%  As a counterpart to the access to the source code and  rights to copy,
%  modify and redistribute granted by the license, users are provided only
%  with a limited warranty  and the software's author,  the holder of the
%  economic rights,  and the successive licensors  have only  limited
%  liability.
%
%  In this respect, the user's attention is drawn to the risks associated
%  with loading,  using,  modifying and/or developing or reproducing the
%  software by the user in light of its specific status of free software,
%  that may mean  that it is complicated to manipulate,  and  that  also
%  therefore means  that it is reserved for developers  and  experienced
%  professionals having in-depth computer knowledge. Users are therefore
%  encouraged to load and test the software's suitability as regards their
%  requirements in conditions enabling the security of their systems and/or
%  data to be ensured and,  more generally, to use and operate it in the
%  same conditions as regards security.
%
%  The fact that you are presently reading this means that you have had
%  knowledge of the CeCILL-B license and that you accept its terms.
%****************************************************************************

%-----------------------------------------------------
%  file : HDF5load.m
%  date : 17/03/2008
%   used to read HDF5 file
%   (developed by Ch. Passeron)
%-----------------------------------------------------
function HDF5load(hdf5_file)

info=hdf5info(hdf5_file);
for i=1:length(info.GroupHierarchy.Datasets)
 name{i}=info.GroupHierarchy.Datasets(i).Name(2:end);
 if (length(info.GroupHierarchy.Datasets(i).Dims)==1)
   data{i}=hdf5read(hdf5_file,name{i})';
 else
   data{i}=hdf5read(hdf5_file,name{i});
 end
 assignin('caller',name{i},data{i}) 
end






