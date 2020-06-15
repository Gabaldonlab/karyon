# cd test
# sh download_data_test_double.sh
# cd ..
echo "Starting Karyon test..."
rm -r ./test/result
python ./bin/2.7/karyon.py -l ./test/DRR040667_1.fastq ./test/DRR040667_2.fastq -d ./test/result --memory_limit 8
