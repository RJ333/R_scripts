#you need: a bio-48 and 49 account (Bernd Schlichting). 48 is accessed via command line, 49 via windows explorer
#a silvaNGS account (if you wanna use this pipeline, this tutorial is for silvaNGS)
#notepad++
#winrar or other unzipping software for .gz-files

#copy raw reads folder (as fastq.gz) to server bio49 --> e.g. to my_folder_on_bio49/miseq/run_1/fastq
#undetermined reads should be excluded, they contain no information
#create more folders in same folder (run_1): joined20, splitreads, cutted 
#if you wanna do this in bash: 
mkdir -p joined20
mkdir -p splitreads
mkdir -p cutted 
#this can cause problems with the rights for editing the folder, our server sees 2 different users when creating folders via windows on bio-49 or on cmd line

#now the data should be where we need it. all necessary programs are already installed (qiime pipeline, cutadapt)
#start mobaxterm
#open ssh-session mit ip-adresse: 10.11.20.48
#login, don't save password
bash			#to enable the bash environment and use commands
#change directory to fastq-folder: 
ll				#"list long" shows all files and folders in your current directory
pwd				#"print working directory" shows the path to your current directory
cd ..			#"change directory" in this case on level above
ll				# with ll we check again what is inside this folder. you can also use 
ls 				#or 
ls -l			#we still need to get one level up
cd ..
ll							#we see a folder named "data"
cd data						#change directory to folder data
cd my_folder_on_bio49		#there we see the folder "my_folder_on_bio49", we could have combined both commands as 
cd data/my_folder_on_bio49  #it now doesn't work anymore, because you are already there
ll							#check what is inside
cd miseq/run_1/fastq		#go to the raw reads

#we found the folder with our reads on bio-48, now we need to join forward and reverse reads. we us the script "multiple_join_paired_ends.py" in QIIME
#flag setzen für minimal read overlap (man page: http://qiime.org/scripts/multiple_join_paired_ends.html)
#parameterfile erstellen: 
vim "dateiname" #opens a text file named "dateiname", press "i" to switch to writing mode and type:
join_paired_end:--min_overlap 20
#press escape to get back to command mode and then :w (for "write") and then :q (for quit)
#now we need to get a level higher (in the run_1 folder)
cd ..
#we start qiime to load the scripts
qiime 
multiple_join_paired_ends.py -i ./fastq/ -o ./joined20 -p ./fastq/dateiname.txt
#when it works you can check via windows in the folder joined20 if new folders are generated there, can take about 2 hours for 96 samples
#in the log_file you will find information on which read was combined with which 
#check files in joined20, you find the joined reads, and the not matching forward and reverse reads. we only use the joined ones
#--> problem: all joined reads in different sample folders have the same name. Here's a small loop that changes the name to the folder's name (incl sample info)
#start this line from the folder joined20
find . -name 'fastqjoin.join.fastq' -exec bash -c 'd="${1%/*}"; mv "$1" "$d/$d-${1##*/}"' - '{}' \; 
#now the joined reads still are in fastq-format, silvaNGS doesn't like that. we have to
#separate quality data from read data with the loop "fastaqual.sh" (must be in correct folder, check with notepad++ or vim or nano) 
#(man page: http://qiime.org/scripts/convert_fastaqual_fastq.html)
bash fastaqual.sh
#this will take a little longer than the joining step, you can check the progress in the splitreads folder
#now we should have all sample files twice, once with .qual ending and once with .fna ending. we need the fna-files and copy those to the folder "cutted"
#in the last preparation step we are going to cut off the primers from the reads, because Illumina only deletes the adapter
#therefore we need our primer sequences
#check sequence files with notepad++ and copy ~25 bp from beginning and end of read
#compare them to your primer sequence (forward as you find it, reverse must be complementary and back->front, sometimes one base is missing or so)
#create primer sequences for the program cutadapt (they need to match the endings as you find them in the reads)
#(man page: https://cutadapt.readthedocs.io/en/stable/)
#we will write a short bash script containing a loop to perform the cutting of the forward primer on all samples
bash cutadapt_primerx_forward.sh
#check some processed files to see the new reads without forward adapter
#all processed files have a slightly changed name, only these should be the input for the second loop --> erase the uncutted files or use different folders
bash cutadapt_primerx_reverse.sh
#if scripts written in notepad++ not work in bash, check the linefeed format (in NP++ "Alle Zeichen anzeigen", dann "bearbeiten -> Zeilenende format --> Unix")
#check some processed files which should now be without any primers
#if some lines are erased completely, silvaNGS will complain during uploading
#it is possible to shorten the primer sequences or search for empty lines and replace with some random bases
#files can now be uploaded to silvaNGS, each upload with max of 500 mb, but multiple uploads parallel are possible
#set your parameters (or use defaults) and request execution
#results: check rarefaction, krona plot, other plots
#download the results archive

#find the csv table results/ssu/tax_breakdown/_otu_breakdown

#this is the final otu table we are going to use for further analyses with R