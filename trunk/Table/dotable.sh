ls -1 *params > input.txt
sed -i 's/.params//g' input.txt
python table.py
