#!/bin/bash
#
#


c=1
c2=1
c3=`ls -d A*/|wc -l`
c4=1

mkdir -p scripts
submission=""
for acc in A*/; do
	
	if [ "$c" -eq 1 ]; then
		submission="scripts/submit_${c2}.sh"
		echo "#!/bin/bash" > $submission
		
	fi
	acc=`basename $acc`;

	echo "tar zcf ${acc}.tar.gz $acc &" >> $submission

	if [ "$c" -eq 32 ] || [ "$c4" -eq "$c3" ];then
		c=1
		echo "wait;" >> $submission
		echo sbatch -J tar_${c2} -t 4:00:00 -A k1616 --output=scripts/slurm-${c2}-%j.err $submission
		c2=$(($c2+1))
	else
		c=$(($c+1))
	fi
	c4=$(($c4+1))
done
