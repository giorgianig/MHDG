# find old values of diffusion
diff_n_old=$(grep -oP 'diff_n =\K.+?(?=! diffusion in the continuity equation)' param.txt)
diff_u_old=$(grep -oP 'diff_u =\K.+?(?=! diffusion in the momentum equation)' param.txt)
diff_e_old=$(grep -oP 'diff_e =\K.+?(?=! perpendicular diffusion for the energy equation)' param.txt)
diff_ee_old=$(grep -oP 'diff_ee =\K.+?(?=! perpendicular diffusion in the electron temperature equation)' param.txt)
diff_nn_old=$(grep -oP 'diff_nn =\K.+?(?=e3 ! diffusion for neutral equation)' param.txt)
diff_vort_old=$(grep -oP 'diff_vort =\K.+?(?=! diffusion in the vorticity equation)' param.txt)

# reset values of diffusion
sed -i "s/diff_n =$diff_n_old/diff_n = $1 /" param.txt
sed -i "s/diff_u =$diff_u_old/diff_u = $1 /" param.txt
sed -i "s/diff_e =$diff_e_old/diff_e = $1 /" param.txt
sed -i "s/diff_ee =$diff_ee_old/diff_ee = $1 /" param.txt
sed -i "s/diff_vort =$diff_vort_old/diff_vort = $1 /" param.txt
sed -i "s/diff_nn =$diff_nn_old/diff_nn = $1/" param.txt
