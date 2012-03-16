for file in ../../../ligands/ZINC/*; do
  vina --config vina.cfg --ligand $file --out out/${file:22:${#file}-22} --log log/${file:22:${#file}-28}.txt
done
