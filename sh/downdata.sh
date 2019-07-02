for i in NAD TAD NSC TSC
do
	echo "------------ Processing $i -------------"
	python py/getcases.py -t $i
done
