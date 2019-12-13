files="
19Bm03.bias.0.40.00
03Am02.dark.1800.36.00
1736342p
07Bm06.flat.z.36.02
03Am05.fringe.z.36.00
2004B.mask.0.36.01
2463796o
1265044o"

obs_ids="
19Bm03.bias.0.40.00
03Am02.dark.1800.36.00
1736342
07Bm06.flat.z.36.02
03Am05.fringe.z.36.00
2004B.mask.0.36.01
2463796
1265044"

for ii in $files; do
  echo $ii
  cadc-data get --fhead --cert $HOME/.ssl/cadcproxy.pem -o ../cfht2caom2/tests/data/$ii.fits.header CFHT $ii.fits.fz || exit $?
done

for ii in $obs_ids; do
  echo $ii
  caom2-repo read --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/ams CFHT $ii > ../cfht2caom2/tests/data/$ii.expected.xml || exit $?
done
