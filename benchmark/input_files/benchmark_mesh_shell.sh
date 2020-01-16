# General commands
# make mesh
python -m salvus_mesher_lite.interface AxiSEMCartesian --basic.model test.bm --basic.period 0.1 --cartesian2Daxisem.x 5 --cartesian2Daxisem.min_z 6366.0 
# run the executable on 4 processors
mpirun -np 4 ./axisem3d 

# Simulation 1 - Cartesian equivalent 
# dir ~/axisem3d_benchmark_sim1
# Sim period 5s
# Cartesian domain radius: 110.97829555
python -m salvus_mesher_lite.interface AxiSEMCartesian --basic.model 3layers.bm --basic.period 5 --cartesian2Daxisem.x 110.97829555 --cartesian2Daxisem.min_z 6321.0 

# Simulation 1 - Spherical equivalent
# dir ~/axisem3d_benchmark_sim2
# Sim period 5s
# Cartesian domain radius: 110.99295622
python -m salvus_mesher_lite.interface AxiSEMCartesian --basic.model 3layers.bm --basic.period 5 --cartesian2Daxisem.x 110.99295622 --cartesian2Daxisem.min_z 6321.0 

# Simulation 1 - Geographic
# dir ~/axisem3d_benchmark_sim3
# Sim period 5s
# Cartesian domain radius: 110.987066515
python -m salvus_mesher_lite.interface AxiSEMCartesian --basic.model 3layers.bm --basic.period 5 --cartesian2Daxisem.x 110.987066515 --cartesian2Daxisem.min_z 6321.0 
