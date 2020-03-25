#!/bin/bash
set -e

current_dir="${PWD}"
remove_part="${SCRAM_ARCH}/cms/cmssw/${CMSSW_VERSION}"
basedir="${CMSSW_RELEASE_BASE/${remove_part}/}"
echo "CMSSW base directory: ${basedir}"


alldirs=()

for DIR in \
  "${ROOTSYS}/include" \
  "${CMSSW_RELEASE_BASE}/src" \
  "$(scram tool tag boost INCLUDE)" \
  "$(scram tool tag tbb INCLUDE)" \
  "$(scram tool tag hepmc HEPMC_BASE)/include" \
  "$(scram tool tag heppdt INCLUDE)" \
  "$(scram tool tag fastjet INCLUDE)" \
  "$(scram tool tag clhep INCLUDE)" \
  "$(scram tool tag python INCLUDE)";
  do
    echo "Adding path: ${DIR/${basedir}/}"
    alldirs+=( "${DIR/${basedir}/}" )
done

pushd "${basedir}"
echo "Creating tarball, may take some time..."
tar -cf "${current_dir}/headers_${CMSSW_VERSION}.tar" "${alldirs[@]}"
echo "headers_${CMSSW_VERSION}.tar was created"
popd
