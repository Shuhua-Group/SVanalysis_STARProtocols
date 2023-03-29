#!/usr/bin/sh

cmdfd=$1
thread=$2

mkfifo $$.fifo
exec 6<>$$.fifo
rm $$.fifo

for((i=0;i<$thread;++i))
do
	echo
done >&6

cat $cmdfd | while read cmd
do
	read -u6 line
	{
		echo $cmd
		$cmd
		echo >&6
	} &
done

wait
echo "All running is ok"

exec 6>&-
