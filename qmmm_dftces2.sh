#!/bin/bash
# DFT-CES2 wrapper script written by Taehwan Jang, Copyright (c) 2024 M-design group @ KAIST
# Global variables
NP=1
QMTYPE="scf" # scf or opt
QMRELAX_HIGH=7 # QM atoms below the specified value are relaxed in unit of Angs
QMMMINISTEP=0 # no of initial QMMM step
QMMMFINSTEP=6 # no of final QMMM step
QMMAXSTEP="5" # QM max trial for "opt"
SUPERCELL=(4 4 1) # (x y z)
SUPERCELLFACTOR=`echo ${SUPERCELL[0]}*${SUPERCELL[1]}*${SUPERCELL[2]} | bc`
DIPOLECORR="yes" # dipole correction for MM charge density
DIPOLEDIR="3" # dipole correction direction x=1 y=2 z=3
DIPOLEPOS="0.80" # dipole correction fractional position along DIPOLEDIR
export DFT_CES2_PATH="/home/jthlol/dft-ces2-bash-test"
LAMMPS="${DFT_CES2_PATH}/MD/lammps/src/lmp_mpi"
QEPW="${DFT_CES2_PATH}/QM/qe/bin/pw.x"
QEPP="${DFT_CES2_PATH}/QM/qe/bin/pp.x"
CUBEADD="${DFT_CES2_PATH}/tools/1-cube/cube_add"
CUBEMULTI="${DFT_CES2_PATH}/tools/1-cube/cube_multi"
BLUR="${DFT_CES2_PATH}/tools/2-blur/blur"
MDDIPOLE="${DFT_CES2_PATH}/tools/3-poisson/dipc2"
LAMMPSIN="base.in.lammps"
QMIN="base.pw.in"
QMIN2="base.pp.in"
LAMMPSRESTART="ini.restart"
MDEQUIL="2000000" # no. of equilibration steps
MDAVERAGE="2000000" # no. of averaging steps
QMLMPTYPE=(3) # i.e. (8 9) type atom type index of QM system defined in LAMMPS data
NSOLVATOM="3000" # i.e. "3000"
MMrepA=(H O)
unit_LAMMPS2QE=0.001686594 # Ry/a.u. to kcal/mol Angs = 592.9107
ATOMNAME=(H Li Be B C N O F Na Mg Al Si P S Cl K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br I W)
ATOMMASS=(1.0 6.9 9.0 10.8 12.0 14.0 15.9 18.9 22.9 24.3 26.9 28.0 30.9 32.0 35.4 39.0 40.0 44.9 47.8 50.9 51.9 54.9 55.8 58.9 58.6 63.5 65.3 69.7 72.6 74.9 78.9 79.9 126.9 183.8)
initialqm=1 #1, when the initial qm has been done.
firstrun=1
cube_coul_QM=(solute.pot.cube) # v_q^A
cube_rho_QM=(solute.rho.cube) # rho^A
cube_QM_rho_hat=(c_rho_hydrogen.cube c_rho_oxygen.cube) # \hat rho^A sigma_b
cube_output_MM=(hydrogen.cube oxygen.cube) # rho_q^B

cube_output_PolarDen=(px.cube py.cube pz.cube) # Polarization density
TIP4Pinvolved=(1 1) # when TIP4P water is involved in the syste, first: 1, second: order of output MM cube for TIP4P oxygen 
if [ ${TIP4Pinvolved[0]} -eq 1 ]; then
	cube_output_TIP4P_O=(rep_TIP4P_oxygen.cube)
fi 
cubeinum=`echo "${#cube_coul_QM[@]}+${#cube_QM_rho_hat[@]}" | bc`
cubeonum=`echo "${#cube_output_MM[@]}+${#cube_output_PolarDen[@]}+${#cube_output_TIP4P_O[@]}" | bc`
cubeio=`echo $cubeinum $cubeonum ${cube_coul_QM[@]} ${cube_QM_rho_hat[@]} ${cube_output_MM[@]} ${cube_output_PolarDen[@]} ${cube_output_TIP4P_O[@]}`
cubeioequil=`echo $cubeinum 0 ${cube_coul_QM[@]} ${cube_QM_rho_hat[@]}`
ATOMrepANAME=(H O N S C Li Na K Rb Cs He Ne Ar Kr Xe F Cl Br I)
ATOMrepA=(1.381 15.56 10.51 3.39 5.15 11.30 24.82 23.51 19.26 11.23 23.78 22.01 11.08 7.617 4.532 2.600 2.694 2.175 1.938)
ATOMalpha=(1.693 5.20 7.25 19.5 11.7 0.193 0.93 5.05 8.32 15.0 1.38 2.67 11.1 16.8 27.2 15.0 30.3 42.8 61.7)
PI=3.141592
temprep=${cube_output_MM}
tempalpha=${cube_outpt_MM}
for ((i=0; i<${#MMrepA[@]}; i++));do
  for ((j=0; j<${#ATOMrepANAME[@]}; j++));do
		if [ "${ATOMrepANAME[$j]}" == "${MMrepA[$i]}" ];then 
      temprep[$i]=${ATOMrepA[$j]}
      tempalpha[$i]=${ATOMalpha[$j]}
    fi
	done
done

Blur_MM (){
	MMpartial=""
	natomtype=`grep 'atom types' $1 | awk '{print $1}'`
	nqmtype=${#QMLMPTYPE}
	nmmtype=$((natomtype-nqmtype))
	line=`grep -n 'Atoms' $1 | cut -f1 -d :`
	for i in $(seq 1 ${nmmtype});do
	temp=`awk -v type=$i -v line=$line 'BEGIN{val=0} {if(NR>line && val==0 && NF==10 && $3==type) {printf("%f\n", $4);val=1;}}' $1`
	MMpartial="$MMpartial${temp} " 
	done
	MMpartialCharge=($MMpartial)

  for ((i=0; i<${#cube_output_MM[@]}; i++));do
		radius=`echo ${PI} ${tempalpha[$i]} | awk '{printf("%f",((2/$1)**0.5*$2/3)**(1/3))}'`
		if [ "$i" == "${TIP4Pinvolved[1]}" ]; then
        $BLUR rep_TIP4P_oxygen.cube ${radius}
			  mv blurred.cube c_rep_TIP4P_oxygen.cube
		else 
			conv=`echo ${MMpartialCharge[$i]} ${temprep[$i]} | awk '{printf("%f",1/$1*$2*2)}'`
			echo $conv
			echo $radius
  	  $CUBEMULTI ${cube_output_MM[$i]} $conv # 1/q*repA/313.7545 
    	$BLUR multiplied.cube ${radius} # atomic radius from isotropic polarizability.
			mv blurred.cube c_rep_${cube_output_MM[$i]}
		fi
	done

        list=$(ls c_rep_*.cube) #convoluted repulsion potential from mm; fixed external potential by mm
        $CUBEADD ${list[@]}
        mv add.cube repA.cube

}
Blur_QM (){
  for ((i=0; i<${#cube_output_MM[@]}; i++));do
			radius=`echo ${PI} ${tempalpha[$i]} | awk '{printf("%f",((2/$1)**0.5*$2/3)**(1/3))}'`
        $BLUR  ${cube_rho_QM} ${radius}
        mv  blurred.cube  ${cube_QM_rho_hat[$i]}
	done
}
post_dipole (){
        dip_parse="${DFT_CES2_PATH}/tools/1-cube/cube_dip_parse" # post-process for polarization density
        addfour="${DFT_CES2_PATH}/tools/1-cube/cube_addfour" 
        $dip_parse ${cube_output_PolarDen[0]} 1          # compute distribution of bound charge along x-direction
        mv postdipole.cube postdipx.cube
        $dip_parse ${cube_output_PolarDen[1]} 2          # compute distribution of bound charge along y-direction
        mv postdipole.cube postdipy.cube
        $dip_parse ${cube_output_PolarDen[2]} 3          # compute distribution of bound charge along z-direction
        mv postdipole.cube postdipz.cube

        # time-averaged charge distribution from MM
        $CUBEADD ${cube_output_MM[@]}
        mv add.cube mobile.cube
        $addfour postdipx.cube postdipy.cube postdipz.cube mobile.cube
        mv add.cube MOBILE_final.cube
        rm add*.cube post*.cube dip?.cube
}

# main loop
for ((qmmmstep=$QMMMINISTEP; qmmmstep<=$QMMMFINSTEP; qmmmstep++));do
  echo "######### Starting $qmmmstep QMMM step #########"$'\n'
  # Make saving directory
  mkdir qm_$qmmmstep mm_$qmmmstep

        ## Preparation of QM
        mpirun -np 1 $LAMMPS -r $LAMMPSRESTART data.temp
  if [ $firstrun -eq 1 ]; then
          # Extract to data from restart
          natoms_qm=`grep nat $QMIN | sed 's/,//' | awk '{print $3}'`
                nqm_extend=`echo ${SUPERCELLFACTOR}*${natoms_qm} | bc `
    # Parsing cell box & atom names
    echo "### Parsing cell box & atom species from data file.."
    Box=(`awk '{if(NF==4 && ($3=="xlo"||$3=="ylo"||$3=="zlo")) print $0}' data.temp`)
    cells=( `printf "%f" ${Box[1]}` `printf "%f" ${Box[5]}` `printf "%f" ${Box[9]}` ) # Cell box sizes 0:x 1:y 2:z
    echo "### Parsed lammps cell box size : ${cells[@]}"$'\n'
    Masses=(`awk 'BEGIN{tag=0;cnt=0}{if($2=="atom"&&$3=="types") natoms=$1; if($1=="Masses") tag=1; if(tag==1 && NF==2 && cnt<natoms) {print $0; cnt++}}' data.temp`)
    for ((i=0; i<${#Masses[@]}; i++));do
      if [ $((i%2)) -eq 1 ]; then
        temp=`echo "${Masses[$i]}*10/1" | bc`  # mass round down to one decimal place
        temp=`echo "$temp/10" | bc -l`
        for ((j=0; j<${#ATOMMASS[@]}; j++));do
          temp2=`echo $temp'=='${ATOMMASS[$j]} | bc -l` # match atomname and atommass
          if [ "$temp2" == "1" ]; then
            atoms="$atoms ${ATOMNAME[$j]}"
          fi
        done
      fi
    done
    atoms=($atoms) # Atomic species by Type in lammps starting from 0..
    echo "### Parsed lammps atomic species : ${atoms[@]}"$'\n' # for all atomic species
  fi # 0th step initial parsing end

  # Parsing solute xyz
  # Parsing dispersion force xyz
  awk -v nsolvatom=$NSOLVATOM '{if(NF==10 && $1>nsolvatom) print $0}' data.temp > data.solute.temp
  forloopend=$((NSOLVATOM+nqm_extend))
  for ((i=${NSOLVATOM}; i<${forloopend}; i++));do
      awk -v var=$i '{if($1==var+1) print $0}' data.solute.temp >> data.solute.temp2
  done
  nodispf=`ls dispf.ave`
  if [ "$nodispf" == "" ]; then
        for ((i=0; i<$nqm_extend; i++));do
                        echo "0 0 0 0 0 0 0 0">> dispf.ave
                done
  fi
  tail -n $nqm_extend dispf.ave | awk -v unit=$unit_LAMMPS2QE '{print $6*unit,$7*unit,$8*unit}' > data.dispforce # sorted id 
  paste data.solute.temp2 data.dispforce > data.solute
  rm data.solute.temp data.solute.temp2
  solute=(`awk -v var="${QMLMPTYPE[*]}" 'BEGIN{split(var,list," ")} {if(NF==13) for(x in list) if($3==list[x]) print $1,$3,$5,$6,$7,$11,$12,$13}' data.solute`) # format: id type x y z bfx bfy bfz, only solute position, x in list: from 1 to the number of QMLMPTYPE
  echo "### Parsed number of atoms in solute : $natoms_qm"

  for ((i=0;i<${#solute[@]};i++));do # i from 0 to 8*natoms; 8 since the element for each atom is id type x y z d_bfx d_bfy d_bfz
    if [ $((i%8)) -eq 1 ]; then
      idx=`echo ${solute[$i]}'-'1 | bc` # since atomic species has allocated from 0; the number is where the atomname is stored.
      solute[$i]=${atoms[$idx]} # swap and allocate the atom type of each solute atom.
    fi
    if [ $((i%8)) -eq 2 -o $((i%8)) -eq 3 -o $((i%8)) -eq 4 -o $((i%8)) -eq 5 -o $((i%8)) -eq 6 -o $((i%8)) -eq 7 ]; then
      solute[$i]=`printf "%lf" ${solute[$i]}`
    fi
  done

  # Solute xyz to qmxyz string
  qmxyz=""
  force=""
  for ((i=0; i<$natoms_qm; i++));do
    qmxyz="$qmxyz${solute[$((8*i+1))]} ${solute[$((8*i+2))]} ${solute[$((8*i+3))]} ${solute[$((8*i+4))]}"$'\n'
                fx=0;   fy=0;   fz=0;
                for ((j=0; j<${SUPERCELLFACTOR}; j++));do
                 fx=`echo ${fx} ${solute[$((8*(i+j*natoms_qm)+5))]} | awk -v super=${SUPERCELLFACTOR} '{print $1+$2/super}' `
                 fy=`echo ${fy} ${solute[$((8*(i+j*natoms_qm)+6))]} | awk -v super=${SUPERCELLFACTOR} '{print $1+$2/super}' `
                 fz=`echo ${fz} ${solute[$((8*(i+j*natoms_qm)+7))]} | awk -v super=${SUPERCELLFACTOR} '{print $1+$2/super}' `
                done
                force="$force${solute[$((8*i+1))]} ${fx} ${fy} ${fz}"$'\n'
  done
  echo "###qmxyz"
  echo "$qmxyz"
  echo "###force"
  echo "$force"

  # Make pw.in based on $QMIN
  awk -v "geo=$qmxyz" -v "dispf=$force" '{if($0=="###qmxyz") {print geo} else if($0=="###dispf") {print dispf} else {print $0} }' $QMIN > pw.in
  if [ $qmmmstep -gt 0 -a $firstrun -eq 0 ]; then 
    sed -i "s/.*\&CONTROL.*/&\ndft_ces = .true./" pw.in
    sed -i "s/.*\&CONTROL.*/&\nrho_ces = '.\/MOBILE_final.cube'/" pw.in
    sed -i "s/.*\&CONTROL.*/&\npauli_rep_ces = '.\/repA.cube'/" pw.in
    sed -i "s/.*\&ELECTRONS.*/&\nstartingwfc = 'file'/" pw.in
    sed -i "s/.*\&ELECTRONS.*/&\nstartingpot = 'file'/" pw.in
  fi
  if [ $qmmmstep -gt 0 -a "$QMTYPE" == "opt" ]; then # relax geometry
		sed -i "s/.*\&CONTROL.*/&\nforc_conv_thr = 1.0D-3/" pw.in
		sed -i "s/.*\&CONTROL.*/&\nnstep = 150/" pw.in
		sed -i "s/.*calculation.*/calculation = 'relax'/" pw.in
		sed -i "s/.*ATOMIC_SPECIES.*/\&IONS\n&/" pw.in
		sed -i "s/.*ATOMIC_SPECIES.*/\/\n&/" pw.in
		sed -i "s/.*\&IONS.*/&\ntrust_radius_max = 0.05D0/" pw.in
		nr=$(grep -n "ATOMIC_POSITIONS" pw.in | cut -d: -f1)
		nrf=$(grep -n "ATOMIC_FORCES" pw.in | cut -d: -f1)
		awk -v line=$nr '{if(NR<=line) {print $0}}' pw.in >tmp
		awk -v line=$nr -v zval=$QMRELAX_HIGH -v nrf=$nrf '{if(NR>line) {if(NF==4 && $4<zval && NR<nrf) {print $0, 0,"", 0, "", 0 } else {print $0}} }' pw.in >>tmp
		mv tmp pw.in
  fi
  firstrun=0
  
  if [ $qmmmstep -gt 0 ]; then
    cp mm_$((qmmmstep-1))/MOBILE_final.cube ./
    cp mm_$((qmmmstep-1))/empty.cube ./
    cp mm_$((qmmmstep-1))/repA.cube ./
  fi
        ## end of Preparation of QM

  # Run QM calculation
  finished=""
  for ((i=0; i<$QMMAXSTEP; i++));do
    if [ $initialqm -eq 0 ]; then
      rm toy
                        mpirun -np $NP $QEPW < pw.in > pw.out
      if [ "$QMTYPE" == "scf" ]; then
        finished=`grep "JOB DONE" pw.out`
        if [ "$finished" != "" ]; then 
          echo "### QM calculation done"$'\n'
          break
        else
          echo "### QM calculation aborted, check output file"$'\n'
          exit
        fi
      elif [ $qmmmstep -eq 0 ]; then
        finished=`grep "JOB DONE" pw.out`
        if [ "$finished" != "" ]; then 
          echo "### QM calculation done"$'\n'
          break
        else
          echo "### QM calculation aborted, check output file"$'\n'
          exit
        fi
      else
        finished=`grep "JOB DONE" pw.out`
        if [ "$finished" != "" ]; then
          echo "### QM calculation done"$'\n'
          break
        else
          echo "### QM calculation aborted, check output file"$'\n'
          exit
        fi
      fi
    fi
  done
		if [ $initialqm -eq 0 ]; then
			solute_nonortho=(`grep -A ${natoms_qm} "ATOMIC_POSITIONS" pw.out | tail -n ${natoms_qm} | awk '{print $1,$2,$3,$4 }'`)
			qmxyz_nonortho=""
			force_nonortho=""
			for ((i=0; i<$natoms_qm; i++));do
				qmxyz_nonortho="$qmxyz_nonortho${solute_nonortho[$((4*i))]} ${solute_nonortho[$((4*i+1))]} ${solute_nonortho[$((4*i+2))]} ${solute_nonortho[$((4*i+3))]}"$'\n'
				force_nonortho="$force_nonortho${solute_nonortho[$((4*i))]} 0 0 0"$'\n'
			done
				awk -v "geo=$qmxyz_nonortho" -v "dispf=$force_nonortho" '{if($0=="###qmxyz") {print geo} else if($0=="###dispf") {print dispf} else {print $0} }' $QMIN > pw.nonortho.in
				sed -i "s/.*\&CONTROL.*/&\ndft_ces = .true./" pw.nonortho.in
				sed -i "s/.*\&CONTROL.*/&\nrho_ces = '.\/MOBILE_final.cube'/" pw.nonortho.in
				sed -i "s/.*\&CONTROL.*/&\npauli_rep_ces = '.\/empty.cube'/" pw.nonortho.in
				mpirun -np $NP $QEPW < pw.nonortho.in > pw.nonortho.out
		fi

    ## Postprocess of QM
  # Generating QM solute potential
  echo "### QM potential post processing.. "$'\n'
  cp $QMIN2 pp.pot.in
  cp $QMIN2 pp.rho.in
  sed -i "s/.*plot_num.*/plot_num = 0/" pp.rho.in
  sed -i "s/.*fileout.*/fileout = 'solute.rho.cube'/" pp.rho.in
 
    if [ $initialqm -eq 0 ]; then
        mpirun -np $NP $QEPP < pp.pot.in > pp.pot.out
        ppdone1=`ls solute.pot.cube`
        mpirun -np $NP $QEPP < pp.rho.in > pp.rho.out
        ppdone2=`ls solute.rho.cube`
        if [ "$ppdone1" == "" ] || [ "$ppdone2" == "" ] ; then
            echo "### QM post processing failed, aborting whole qmmm loop"$'\n'
        else
        Blur_QM solute.rho.cube
        $CUBEADD v_saw.cube solute.pot.cube
        $CUBEMULTI solute.pot.cube 0
        mv multiplied.cube empty.cube
        mv add.cube solute.pot.cube
        cp pw.in pw.out pw.nonortho.out pp.pot.out pp.rho.out solute.rho.cube solute.pot.cube qm_$qmmmstep
        cp -r solute qm_$qmmmstep 
        fi
    else
        rhoexst=`ls solute.rho.cube`
        emptyexst=`ls empty.cube`
        if [ "$rhoexst" != "" ] ; then
            Blur_QM solute.rho.cube
        else
            cp qm_$qmmmstep/*.cube .
            Blur_QM solute.rho.cube
        fi
        if [ "$emptyexst" == "" ] ; then
            $CUBEMULTI solute.rho.cube 0
            mv multiplied.cube empty.cube
        fi
    fi
    initialqm=0
    echo "qmmstep is $qmmmstep"
    # Updating QM optimized geometry
    cp qm_$qmmmstep/pw.out .
    if [ $qmmmstep -gt 0 -a "$QMTYPE" == "opt" ]; then
        solutefinal=(`grep -A ${natoms_qm} "ATOMIC_POSITIONS" pw.out | tail -n ${natoms_qm} | awk '{print $1,$2,$3,$4 }'`)
        cnt=0
        for ((i=0; i<${#solute[@]}; i++));do
            if [ $((i%8)) != 0 -a $((i%8)) != 1 -a $((i%8)) != 5 -a $((i%8)) != 6 ]; then
                solute[$i]=${solutefinal[$((i-4*cnt-1))]}
            fi
            if [ $((i%8)) == 7 ]; then let cnt=cnt+1; fi
        done
        echo "### Updating QM optimized geometry.. "$'\n'
    fi
    ## end of Postprocess of QM

    ## Preparation of MM
    # Modifying LAMMPS INPUT : restart, geometry
    cp $LAMMPSIN in.lammps
    if [ $qmmmstep -gt 0 ]; then
			for ((i=0; i<$natoms_qm; i++));do
				cnt=0
				for ((j=1; j<=${SUPERCELL[0]}; j++));do
					for ((k=1; k<=${SUPERCELL[1]}; k++));do
						for ((l=1; l<=${SUPERCELL[2]}; l++));do
							index=`echo "${solute[$((i*8))]} + $cnt*$natoms_qm" | bc`
							xpos=`echo ${solute[$((i*8+2))]} ${cells[0]} $((j-1)) ${SUPERCELL[0]} | awk '{print $1+$2*$3/$4}'`
							ypos=`echo ${solute[$((i*8+3))]} ${cells[1]} $((k-1)) ${SUPERCELL[1]} | awk '{print $1+$2*$3/$4}'`
							zpos=`echo ${solute[$((i*8+4))]} ${cells[2]} $((l-1)) ${SUPERCELL[2]} | awk '{print $1+$2*$3/$4}'`
							cnt=$((cnt+1))
							sed -i "s/.*xyz.*/&\nset atom ${index} x ${xpos} y ${ypos} z ${zpos}/" in.lammps
						done
					done
				done
			done
    fi
    ## end of Preparation of MM

    # RUN LAMMPS equilibration step
    mpirun -np 1 $LAMMPS -r $LAMMPSRESTART data.md
    sed -i "s/.*read_data.*/read_data data.md/" in.lammps
    echo "### Running LAMMPS(emxext) for equilibration $qmmmstep QMMM iterations"$'\n'
    sed -i "s/.*run.*/run\t\t$MDEQUIL/" in.lammps
    if [ $qmmmstep -eq 0 ]; then
        LAMMPSRESTARTtime=0
    else
        LAMMPSRESTARTtime=`ls -lrt *.restart | tail -n 1 | awk '{print $9}' | cut -f2 -d '.'`
    fi
  sed -i "s/.*reset_timestep.*/reset_timestep ${LAMMPSRESTARTtime}/" in.lammps
  cp in.lammps in.lammps.equil 
  finaltime=`echo ${LAMMPSRESTARTtime} ${MDEQUIL} | awk '{print $1+$2}'`
  sed -i "s/.*c_fdisp.*/fix\t\tshowf all ave\/atom 1 $MDEQUIL $finaltime c_fdisp[*]/" in.lammps.equil # MDEQUIL is a divisor of MDAVERAGE.
  sed -i "s/.*dispf.ave.*/dump\t\t11 SOLUTE custom $finaltime dispf.ave id type xu yu zu f_showf[*]/" in.lammps.equil
sed -i "s/.*#CUBEPOSITION.*/&\ngrid\t\t ${cubeioequil} /" in.lammps.equil
  mpirun -np $NP $LAMMPS -in in.lammps.equil > lammps.equil.out

  # RUN LAMMPS averaging step
  LAMMPSRESTART=`ls -lrt *.restart | tail -n 1 | awk '{print $9}'`
  LAMMPSRESTARTtime=`ls -lrt *.restart | tail -n 1 | awk '{print $9}' | cut -f2 -d '.'`
  sed -i "s/.*reset_timestep.*/reset_timestep ${LAMMPSRESTARTtime}/" in.lammps
sed -i "s/.*#CUBEPOSITION.*/&\ngrid\t\t ${cubeio} /" in.lammps
  mpirun -np 1 $LAMMPS -r $LAMMPSRESTART data.md
  finaltime=`echo ${LAMMPSRESTARTtime} ${MDAVERAGE} | awk '{print $1+$2}'`
  sed -i "s/.*read_data.*/read_data data.md/" in.lammps
  sed -i "s/.*c_fdisp.*/fix\t\tshowf all ave\/atom 1 $MDAVERAGE $finaltime c_fdisp[*]/" in.lammps # MDEQUIL is a divisor of MDAVERAGE.
  sed -i "s/.*dispf.ave.*/dump\t\t11 SOLUTE custom $finaltime dispf.ave id type xu yu zu f_showf[*]/" in.lammps # MDEQUIL is a divisor of MDAVERAGE.
  sed -i "s/.*run.*/run\t\t$MDAVERAGE/" in.lammps
  echo "### Running LAMMPS(emdext) for averaging solvent charge density $qmmmstep QMMM iterations"$'\n'
  mpirun -np $NP $LAMMPS -in in.lammps > lammps.average.out
  LAMMPSRESTART=`ls -lrt *.restart | tail -n 1 | awk '{print $9}'`
  rm -f log.lammps
  post_dipole
	Blur_MM
  if [ "$DIPOLECORR" == "yes" ]; then
    $MDDIPOLE MOBILE_final.cube $DIPOLEDIR $DIPOLEPOS
  fi
  cp dispf.ave repA.cube in.lammps* empty.cube MOBILE_final.cube *.lammpstrj $LAMMPSRESTART lammps.*.out mm_$qmmmstep

done # qmmm loop



