cmake_minimum_required(VERSION 2.8)

project(TRIPHIC)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_utils/CMakeModules/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -W -Wall  -Wshadow -Wno-error=unused-function -fomit-frame-pointer -ffast-math -fvisibility=hidden -fvisibility-inlines-hidden -fPIC  -Werror " )

if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^3.1" )
  message( "Hopefully, the Ubuntu 14.04 machine, with a modern compiler." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()
if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6." )
  message( "It looks like a trusty old Centos 6 machine with an ancient compiler." )
endif()

#############################################################################
## required packages
#############################################################################

find_package(OEToolkits COMPONENTS oeszybki oeff oeieff oeamber oequacpac
  oeam1bcc oeam1 oeshape
  oespicoli oesheffield oemmff oemolpotential oeopt oebio oezap
  oegrid oefizzchem oedepict oechem oesystem oeplatform oeiupac oeshape
  oeomega2)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options regex filesystem system)
message( "Boost include dir : ${Boost_INCLUDE_DIRS}" )
message( "Boost version : ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}" )

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_FILE_OFFSET_BITS=64 -ffast-math -Wall -mssse3 -msse2" )
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg" )

#set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg" )
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -D_FILE_OFFSET_BITS=64 -fomit-frame-pointer -ffast-math -Wall -mssse3 -msse2" )

include_directories( SYSTEM ${OEToolkits_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})

#############################################################################
## build triphic
#############################################################################


# differentiate between DEBUG and RELEASE exes
set(EXECUTABLE_OUTPUT_PATH ${TRIPHIC_SOURCE_DIR}/../exe_${CMAKE_BUILD_TYPE})

set(TRIPHIC_SRCS
  triphic.cc
  triphic_subs.cc
  find_cliques.cc
  get_given_oeconf.cc 
  grid_volumes.cc
  make_pphore_sites.cc
  oemol_conf_name_fns.cc
  overlay_oemolbase.cc
  score_and_store_cliques.cc
  send_messages_to_tp_slaves_mpi.cc
  send_pphore_sites_by_mpi.cc
  step_oemolstream.cc
  AtomSphere.cc
  BasePPhoreSite.cc
  OverlayScore.cc
  ParseSMARTSXML.cc
  ShapeTanimoto.cc
  SinglePPhoreSite.cc
  TriphicSettings.cc
  VolumeGrid.cc
  VolumeGridMPI.cc
  build_sites_from_moe_features.cc
  build_vol_from_moe_volumesphere.cc
  read_moe_ph4_file.cc
  read_site_from_file.cc)

set(TRIPHIC_DACLIB_SRCS Combinator.cc
  DiamondOverlay.cc
  PharmPoint.cc
  OverlayTrans.cc
  angle_about_axis.cc
  atomic_num_to_rad.cc 
  bron_kerbosch_mk4.c
  build_time.cc
  centre_coords.cc
  create_cansmi.cc
  eigen_solve_4b4f.cc
  mpi_string_subs.cc
  send_environment_var_to_mpi_slave.cc
  split_filename.cc
  read_smarts_file.cc
  rotate_about_axis.cc)

set(TRIPHIC_INCS gtpl_defs.H
  AtomSphere.H
  BasePPhoreSite.H
  MOEPointsDefs.H
  OverlayScore.H
  ParseSMARTSXML.H
  ShapeTanimoto.H
  SinglePPhoreSite.H
  TriphicSettings.H
  VolumeGrid.H
  VolumeGridMPI.H)

set(TRIPHIC_DACLIB_INCS stddefs.H
  Combinator.H
  DefaultPointsDefs.H
  DiamondOverlay.H
  FileExceptions.H
  PharmPoint.H
  OverlayTrans.H
  crash.H
  points_on_sphere.H)

# make the executable line for triphic - include the include files so that
# QtCreator picks them up and puts them in the project
add_executable(triphic ${TRIPHIC_SRCS} ${TRIPHIC_DACLIB_SRCS}
  ${TRIPHIC_INCS} ${TRIPHIC_DACLIB_INCS})
target_link_libraries(triphic ${LIBS} ${TRIPHIC_SOURCE_DIR}/libtautenum.a ${OEToolkits_LIBRARIES}
  ${Boost_LIBRARIES} ${MPI_LIBRARIES} z pthread)

#############################################################################
## build loob
#############################################################################

set(LOOB_SRCS loob.cc
  loob_subs.cc
  make_pphore_sites.cc
  oemol_conf_name_fns.cc
  BasePPhoreSite.cc
  LoobSettings.cc
  SinglePPhoreSite.cc)

set(LOOB_DACLIB_SRCS angle_about_axis.cc
  atomic_num_to_rad.cc
  build_time.cc
  centre_coords.cc
  eigen_solve_4b4f.cc
  read_smarts_file.cc
  rotate_about_axis.cc
  DiamondOverlay.cc
  OverlayTrans.cc
  ParseSMARTSXML.cc
  PharmPoint.cc)

set(LOOB_DACLIB_INCS stddefs.H
  BasePPhoreSite.H
  DefaultPointsDefs.H
  FileExceptions.H
  ParseSMARTSXML.H
  PharmPoint.H
  LoobSettings.H
  SinglePPhoreSite.H
  SMARTSExceptions.H)

# make the executable line for loob - include the include files so that
# QtCreator picks them up and puts them in the project
add_executable(loob ${LOOB_SRCS} ${LOOB_DACLIB_SRCS}
 ${LOOB_INCS} ${LOOB_DACLIB_INCS})
target_link_libraries(loob ${LIBS} ${OEToolkits_LIBRARIES}
 ${Boost_LIBRARIES} z pthread)

#############################################################################
## build plurality
#############################################################################

set(PLURALITY_SRCS plurality.cc
  plurality_subs.cc
  AtomSphere.cc
  BasePPhoreSite.cc
  PluralityHit.cc
  PluralitySettings.cc
  PPhoreQuery.cc
  SinglePPhoreSite.cc
  VolumeGrid.cc
  VolumeGridMPI.cc
  get_given_oeconf.cc
  grid_volumes.cc
  make_pphore_sites.cc
  oemol_conf_name_fns.cc
  overlay_oemolbase.cc
  read_site_from_file.cc
  send_messages_to_tp_slaves_mpi.cc
  step_oemolstream.cc)

set(PLURALITY_DACLIB_SRCS DiamondOverlay.cc
  ParseSMARTSXML.cc
  PharmPoint.cc
  OverlayTrans.cc
  angle.cc
  angle_about_axis.cc
  atomic_num_to_rad.cc
  build_time.cc
  centre_coords.cc
  eigen_solve_4b4f.cc
  mpi_string_subs.cc
  send_environment_var_to_mpi_slave.cc
  send_oemols_via_mpi.cc
  split_filename.cc
  read_smarts_file.cc
  rotate_about_axis.cc
  torsion.cc)

set(PLURALITY_INCS AtomSphere.H
  BasePPhoreSite.H
  PluralityHit.H
  PluralitySettings.H
  PPhoreQuery.H
  SinglePPhoreSite.H
  VolumeGrid.H
  VolumeGridMPI.H)

set(PLURALITY_DACLIB_INCS DiamondOverlay.H
  PharmPoint.H
  OverlayTrans.H)

# make the executable line for plurality - include the include files so that
# QtCreator picks them up and puts them in the project
add_executable(plurality ${PLURALITY_SRCS} ${PLURALITY_DACLIB_SRCS}
 ${PLURALITY_INCS} ${PLURALITY_DACLIB_INCS})
target_link_libraries(plurality ${LIBS} ${OEToolkits_LIBRARIES}
 ${Boost_LIBRARIES} ${MPI_LIBRARIES} z pthread)

