# cd test
# sh test.sh
# cd ..
echo "Starting Karyon test..."
rm -r ./test/result
python ./bin/2.7/karyon.py -l ./test/DRR040667.fastq -d ./test/result