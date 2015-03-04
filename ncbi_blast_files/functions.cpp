void Write_GPU_information_file_preliminary(FILE* GPU_information_file,
                                            const Int4 num_blocksx,
                                            const Int4 num_threadsx,
                                            const char percentage,
                                            const Int4 num_volumes)
{
    fprintf(GPU_information_file,"#Number of blocksx\n");
    fprintf(GPU_information_file,"%d\n", num_blocksx );
    fprintf(GPU_information_file,"#Number of threadsx\n");
    fprintf(GPU_information_file,"%d\n", num_threadsx);
    fprintf(GPU_information_file,"#Percentage\n");
    fprintf(GPU_information_file,"%d\n", percentage);
    fprintf(GPU_information_file,"#Number of volumes\n");
    fprintf(GPU_information_file,"%d\n", num_volumes );

}
void Write_GPU_information_file(FILE* GPU_information_file,
                                const char* GPU_Database_name,
                                const unsigned long GPU_BLASTP_database_bytes,
                                const Int4 volume_length,
                                const Int4 num_sequences_to,
                                const Int4 store_limit,
                                const Int4 break_limit,
                                const Int4 Group_number,
                                const Int4* Sequence_length_vector){

    fprintf(GPU_information_file,"#Volume name\n");
    fprintf(GPU_information_file,"%s\n", GPU_Database_name);
    fprintf(GPU_information_file,"#Volume size in bytes\n");
    fprintf(GPU_information_file,"%ld\n", GPU_BLASTP_database_bytes );
    fprintf(GPU_information_file,"#Volume sequences\n");
    fprintf(GPU_information_file,"%d\n", volume_length );
    fprintf(GPU_information_file,"#Number of sequences in GPU volume\n");
    fprintf(GPU_information_file,"%d\n", num_sequences_to );
    fprintf(GPU_information_file,"#Store limit\n");
    fprintf(GPU_information_file,"%d\n", store_limit );
    fprintf(GPU_information_file,"#Break limit\n");
    fprintf(GPU_information_file,"%d\n", break_limit );
    fprintf(GPU_information_file,"#Group number\n");
    fprintf(GPU_information_file,"%d\n", Group_number );
    fprintf(GPU_information_file,"#Length of each group\n");
    Int4 i = 0;
    for(i = 0; i < Group_number; ++i)
        fprintf(GPU_information_file,"%d\n", Sequence_length_vector[i] );

}

void Read_GPU_information_file_preliminary(FILE** GPU_information_file,
                                           const BlastSeqSrc* seq_src,
                                           Int4* num_blocksx, Int4* num_threadsx,
                                           Int4* percentage, Int4* num_volumes)
{
    char* Database_name = BlastSeqSrcGetName(seq_src);
    char* GPU_information_name = (char*) calloc(strlen(Database_name) + strlen(".gpuinfo") + 1,sizeof(char));
    strcat(GPU_information_name, Database_name);
    strcat(GPU_information_name, ".gpuinfo");

    *GPU_information_file = fopen( GPU_information_name, "rb");
    if( NULL == (*GPU_information_file) ) {
        printf("GPU information file cannot be opened for reading. Exiting...\n");
        exit(1);
    }

    char *line_read = (char *) malloc( sizeof(char) );
    size_t nbytes = 1;
    Int4 bytes_read = getline(&line_read, &nbytes, *GPU_information_file);//Read first line which is a comment
    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);//Read number of blocksx
    *num_blocksx = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);//Read number of threads
    *num_threadsx = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);//Read percentage
    *percentage = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);//Read number of volumes
    *num_volumes = atoi(line_read);

    free(line_read);
    free(GPU_information_name);

}

void Read_GPU_information_file(const FILE* GPU_information_file,
                               char** GPU_volume_name,
                               unsigned long* GPU_BLASTP_database_bytes,
                               Int4* volume_length,
                               Int4* num_sequences_to,
                               Int4* store_limit,
                               Int4* break_limit,
                               Int4* Group_number,
                               Int4** Sequence_length_vector,
			       const Int4 cpu_threads){

    char *line_read = (char *) malloc( sizeof(char) );
    size_t nbytes = 1;
    Int4 bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read first line which is a comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the volume name
    char* GPU_volume_name_temp = (char*)realloc(*GPU_volume_name, bytes_read*sizeof(char) );
    if( NULL != GPU_volume_name_temp ){
        *GPU_volume_name = GPU_volume_name_temp;
        memset(*GPU_volume_name, 0, bytes_read*sizeof(char));
        memcpy(*GPU_volume_name, line_read, bytes_read - 1);
    }
    else {
        printf("Error (re)allocating memory for the GPU information file\n");
        exit(1);
    }

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the size of the volume in bytes
    *GPU_BLASTP_database_bytes = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the number of sequences in this volume
    *volume_length = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the number of sequences processed by the GPU in this volume
    *num_sequences_to = atoi(line_read);

    switch ( cpu_threads ){

    case 1:
      *num_sequences_to = ((*volume_length) * 95) / 100;
      break;
    case 2:
      *num_sequences_to = ((*volume_length) * 85) / 100;
      break;
    case 3:
      *num_sequences_to = ((*volume_length) * 80) / 100;
      break;
    case 4:
      *num_sequences_to = ((*volume_length) * 76) / 100;
      break;
    case 5:
      *num_sequences_to = ((*volume_length) * 73) / 100;
      break;
    case 6:
      *num_sequences_to = ((*volume_length) * 70) / 100;
      break;
    default:
      *num_sequences_to = ((*volume_length) * 70) / 100;
      break;

    }

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the ID of the last sequences processed by the GPU in this volume
    *store_limit = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the ID of the last sequences of this volume
    *break_limit = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the number of groups in this volume
    *Group_number = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    Int4* Sequence_length_vector_temp = (Int4*) realloc(*Sequence_length_vector, (*Group_number)*sizeof(Int4));
    if( NULL != Sequence_length_vector_temp )
        *Sequence_length_vector = Sequence_length_vector_temp;
    else {
        printf("Error (re)allocating memory for the GPU information file\n");
        exit(1);
    }
    Int4 i = 0;
    for( i = 0; i < *Group_number; ++i){
        bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the length of each group of this volume
        (*Sequence_length_vector)[i] = atoi(line_read);
    }

    free(line_read);

}

int ExtractVolumeSequences(unsigned char * buffer){
    //Extracting information from the index file (see index_file.txt for the format information)

    int version =
        buffer[0 + 0]*(1<<24) +
        buffer[0 + 1]*(1<<16) +
        buffer[0 + 2]*(1<<8) +
        buffer[0 + 3];
    int type =
        buffer[4 + 0]*(1<<24) +
        buffer[4 + 1]*(1<<16) +
        buffer[4 + 2]*(1<<8) +
        buffer[4 + 3];
    if( type != 1 ){
        fprintf(stderr,"Error index file is not a protein index file\n");
        exit(1);
    }
    int title_bytes =
        buffer[8 + 0]*(1<<24) +
        buffer[8 + 1]*(1<<16) +
        buffer[8 + 2]*(1<<8) +
        buffer[8 + 3];
    int date_bytes =
        buffer[8 + 4 + title_bytes + 0 ]*(1<<24) +
        buffer[8 + 4 + title_bytes + 1 ]*(1<<16) +
        buffer[8 + 4 + title_bytes + 2 ]*(1<<8) +
        buffer[8 + 4 + title_bytes + 3 ];
    int sequences =
        buffer[8 + 4 + 4 + title_bytes + date_bytes + 0 ]*(1<<24) +
        buffer[8 + 4 + 4 + title_bytes + date_bytes + 1 ]*(1<<16) +
        buffer[8 + 4 + 4 + title_bytes + date_bytes + 2 ]*(1<<8) +
        buffer[8 + 4 + 4 + title_bytes + date_bytes + 3 ];

    return sequences;
}

void ReadGPU_BLASTP_volinfo(const BlastSeqSrc* seq_src, Int4 *num_volumes, Int4** volume_lengths)
{
    if( GPU_VERBOSE )
        printf("Reading volume information\n");

    Int4* volume_lengths_temp;
    *volume_lengths = NULL;

    //deep copy the database path to a temporary string because "dirname()" might modify the database path
    char* Database_name = BlastSeqSrcGetName(seq_src);
    int Database_name_bytes = strlen( Database_name ) + 1;
    char* Database_copy = (char*) calloc( Database_name_bytes, sizeof(char) );
    memcpy( Database_copy, Database_name, Database_name_bytes );

    char *aliasfile_name = (char*)calloc(strlen(Database_copy) + strlen(".pal") + 1, sizeof(char));
    strcat(aliasfile_name, Database_copy);
    strcat(aliasfile_name, ".pal");
    FILE* aliasfile = fopen(aliasfile_name, "r" );

    if( NULL == aliasfile ) {//if the the protein alias file (.pal) file doesn't exist then

        *num_volumes = 1; //there is only one volume

        char *indexfile_name = (char*)calloc(strlen(Database_copy) + strlen(".pin") + 1, sizeof(char));
        strcat(indexfile_name, Database_copy);
        strcat(indexfile_name, ".pin");
        FILE* indexfile = fopen(indexfile_name, "r" );//read index file

        if( NULL == indexfile) {
            fprintf(stderr, "Cannot open index file %s\n", indexfile_name);
            exit(1);
        } else {
            size_t result = 0;
            unsigned char* buffer = (char*)malloc(10000*sizeof(char));
            result = fread(buffer, sizeof(char), 10000, indexfile);//read the first 10000 bytes of the index file

            volume_lengths_temp = (Int4*) realloc( *volume_lengths, 1*sizeof(Int4) );
            if( NULL != volume_lengths_temp )
                *volume_lengths = volume_lengths_temp;
            else {
                printf("Error (re)allocating during index file reading\n");
                exit(1);
            }

            (*volume_lengths)[0] = ExtractVolumeSequences(buffer);

            free(buffer);
        }

        free(indexfile_name);
        fclose(indexfile);

    } else { //else read the lines line by line to find the "DBLIST" line

        char* line_read = (char *) malloc( sizeof(char) );
        size_t nbytes = 1;
        Int4 bytes_read = getline(&line_read, &nbytes, aliasfile);//Read firs line

        char * pch;
        *num_volumes = 0;
        int k = 0;
        for(k = 0; k < 1000; ++k){//read the first 1000 lines of the alias file

            if( feof(aliasfile) ) //at end of file break
                break;

            pch = strtok (line_read," "); //extract the first word of the line
            int diff = strncmp(pch, "DBLIST", 6); //compare the word with "DBLIST"
            if(diff){ //if it is not "DBLIST" read the next line and continue

                bytes_read = getline(&line_read, &nbytes, aliasfile);//Read firs line which is a comment
                continue;

            } else { //if it is "DBLIST" read the rest of the volume names

                pch = strtok (NULL, " "); //extract the first volume name
                while (pch != NULL){

                    //check if the last character of pch is '\n' which means this is the last index file
                    if( pch[strlen(pch) - 1 ] == '\n')
                        pch[strlen(pch) - 1 ] = '\0'; //replace with the null character

                    char *indexfile_name = (char*)calloc(strlen(pch) + strlen(".pin") + 1, sizeof(char));
                    strcat(indexfile_name, pch);
                    strcat(indexfile_name, ".pin");
                    FILE* indexfile = fopen(indexfile_name, "r" );//read index file

                    if( NULL == indexfile){

                        fprintf(stderr, "Cannot open the index file %s\n", indexfile_name);
                        exit(1);

                    } else {
                        size_t result = 0;
                        unsigned char* buffer = (char*)malloc(10000*sizeof(char));
                        result = fread(buffer, sizeof(char), 1000, indexfile);//read the first 1000 bytes of the index file


                        volume_lengths_temp = (Int4*) realloc( *volume_lengths, ((*num_volumes)+1)*sizeof(Int4) );
                        if( NULL != volume_lengths_temp )
                            *volume_lengths = volume_lengths_temp;
                        else {
                            printf("Error (re)allocating during index file reading\n");
                            exit(1);
                        }

                        (*volume_lengths)[*num_volumes] = ExtractVolumeSequences(buffer);

                        printf("%s: %d\n",indexfile_name, (*volume_lengths)[*num_volumes]);
                        free(buffer);

                    }

                    (*num_volumes)++;

                    free(indexfile_name);
                    fclose(indexfile);
                    pch = strtok (NULL, " "); //extract the next volume name

                }


                break; //after done with reading all volume names break
            }

        }

        free(line_read);
    }

    if( NULL != aliasfile )
        fclose(aliasfile);

    free(aliasfile_name);
    free(Database_copy);

}


/** Creates the database used by the GPU_BLASTP
 @param Sequence_length_vector vector that holds the length of each sequence group [in]
 @param Group_number           the number of groups that the sequence database is split [in]
 @param Group_size             size of each group of sequences [in]
 @param stride                 that is used in createing the GPU database [in]
 @param seq_arc                data structure that holds information about the sequence database [in]
 @param num_sequence_to        last sequence of the database that is stored in the GPU database [in]
*/
void CreateGPU_BLASTP_database(BlastSeqSrc* seq_src,
                               const BlastGPUOptions* gpu_options)
{


  Int4 i = 0, j = 0;
  Int4 stride = 4;
  Int4 Group_size = (gpu_options->num_blocksx) * (gpu_options->num_threadsx);
  char percentage = PERCENTAGE;

  BlastSeqSrcIterator* itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
  BlastSeqSrcGetSeqArg seq_arg;
  memset((void*) &seq_arg, 0, sizeof(seq_arg));

  /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
     sequences are retrieved in ncbistdaa/ncbi2na encodings respectively. */
  seq_arg.encoding = eBlastEncodingProtein;

  char* Database_name = BlastSeqSrcGetName(seq_src);
  char* GPU_Database_name = NULL;
  char* GPU_Database_name_temp;

    char* GPU_BLASTP_database = NULL;
    char* GPU_BLASTP_database_temp;
    Int4* Sequence_length_vector = NULL;
    Int4* Sequence_length_vector_temp;
    char* digits = (char*) calloc(10, sizeof(char));//10 digits are enough to cover up to 999999999 volumes

    char* GPU_information_name = (char*) calloc(strlen(Database_name) + strlen(".gpuinfo") + 1,sizeof(char));
    strcat(GPU_information_name, Database_name);
    strcat(GPU_information_name, ".gpuinfo");

    FILE* GPU_information_file = fopen( GPU_information_name, "wb");
    if ( NULL == GPU_information_file )
        printf("GPU information file cannot be opened for writing\n");

    Int4* volume_lengths = NULL;
    Int4 num_volumes = 0;
    ReadGPU_BLASTP_volinfo(seq_src, &num_volumes, &volume_lengths);
    printf("num_volumes = %d\n",num_volumes);
    Write_GPU_information_file_preliminary(GPU_information_file, gpu_options->num_blocksx,
                                           gpu_options->num_threadsx, percentage, num_volumes);

    Int4 volume = 0, break_limit = 0;
    for(volume = 0; volume < num_volumes; ++volume)
    {
        Int4 num_sequences_to = (percentage * volume_lengths[volume]) / 100;
        Int4 store_limit = break_limit + num_sequences_to;
        Int4 Group_number = (num_sequences_to / Group_size ) + (( (num_sequences_to % Group_size ) == 0)?0:1);

        Sequence_length_vector_temp = (Int4*) realloc( Sequence_length_vector, Group_number*sizeof(Int4) );
        if( NULL != Sequence_length_vector_temp )
            Sequence_length_vector = Sequence_length_vector_temp;
        else {
            printf("Error (re)allocating memory for the GPU database creation\n");
            exit(1);
        }

        Int4 vector_counter = 0;
        unsigned long GPU_BLASTP_database_bytes = 0;
        Int4 oid = break_limit + Group_size - 1; //initialized oid to the last sequence of the first group
        if( oid >= store_limit )
            oid = store_limit - 1;

        //Initializes the Sequence_length_vector with the length of the longest sequence in each group
        //The length is converted to a multiple of WIDTH_MULTIPLE (defined in gpu_cpu_common.h)
	int previous_sequence_length = 0;
        do {

            Int4 sequence_length = BlastSeqSrcGetSeqLen( seq_src, (void*) &oid );
 	    if( previous_sequence_length > sequence_length ){
 	      fprintf(stderr,"ERROR: the input database is not sorted\n");
 	      fprintf(stderr,"Sort the database by using the option \"-sort_volumes\" with \"makeblastdb.\"\n");
 	      fprintf(stderr,"Exiting...\n");
	      exit(1);
 	    }
	    previous_sequence_length = sequence_length;

            Round2Multiple( &sequence_length );
            Sequence_length_vector[ vector_counter ] = sequence_length;
            if( GPU_VERBOSE )
                printf("Group %d starts at sequence %d and has length = %d\n",
                       vector_counter, oid, Sequence_length_vector[vector_counter]);

            ++vector_counter;

            if( oid == store_limit - 1 )
                break;

            oid += Group_size;

        } while( oid < store_limit - 1 );

        //If the last iteration of the previous do-while didn't cover all sequences (i.e. the Group_size is not a multiple of num_sequences)
        //the remaining sequences form the last group, and the last element of Sequence_length_vector is initialized
        if( oid != store_limit - 1 ){
            oid = store_limit - 1;
            Int4 sequence_length = BlastSeqSrcGetSeqLen( seq_src, (void*) &oid );
	    /* check if the volume is sorted */
	    if( previous_sequence_length > sequence_length ){
 	      fprintf(stderr,"ERROR: the input database is not sorted\n");
 	      fprintf(stderr,"Sort the database by using the option \"-sort_volumes\" with \"makeblastdb.\"\n");
 	      fprintf(stderr,"Exiting...\n");
	      exit(1);
 	    }
            Round2Multiple( &sequence_length );
            Sequence_length_vector[ vector_counter ] = sequence_length;
            if( GPU_VERBOSE )
                printf("Group %d starts at sequence %d and has length = %d\n",
                       vector_counter, oid, Sequence_length_vector[vector_counter]);
        }

        //Calculate the total size of the GPU database
        for(i = 0; i < Group_number; ++i)
            GPU_BLASTP_database_bytes += Sequence_length_vector[i]*Group_size*sizeof(char);

        //This helps to check if there is a very large sequence that can cause the GPU memory to explode
        if( GPU_BLASTP_database_bytes > UINT4_MAX )
        {
            printf("Memory requirements for GPU database larger than 4 GB. Exiting...\n");
            exit(0);
        }

        GPU_BLASTP_database_temp = (char*) realloc( GPU_BLASTP_database, GPU_BLASTP_database_bytes*sizeof(char) );
        if( NULL != GPU_BLASTP_database_temp )
            GPU_BLASTP_database = GPU_BLASTP_database_temp;
        else {
            printf("Error (re)allocating memory for the GPU database creation\n");
            exit(1);
        }

        memset(GPU_BLASTP_database, 0, GPU_BLASTP_database_bytes*sizeof(char));

        int length = 0;
        int group_member = 0;
        int counter = 0, sequence = 0, AminoAcids = 0, base = 0;
        break_limit += volume_lengths[volume];
        /* iterate over all subject sequences of this volume */
	previous_sequence_length = 0;
        while ( ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) ) {

            if (seq_arg.oid == BLAST_SEQSRC_ERROR)
                break;
            if (BlastSeqSrcGetSequence(seq_src,  &seq_arg) < 0)
                continue;


            if( seq_arg.oid < store_limit ){

                if( (counter % Group_size) == 0 ){
                    base += length*Group_size;
                    length = Sequence_length_vector[ group_member ];
                    group_member++;
                    sequence = 0;
                }

                AminoAcids = 0;
                int sequence_length = BlastSeqSrcGetSeqLen( seq_src, (void*) &seq_arg.oid );
		if( previous_sequence_length > sequence_length ){
		  fprintf(stderr,"ERROR: the input database is not sorted\n");
		  fprintf(stderr,"Sort the database by using the option \"-sort_volumes\" with \"makeblastdb.\"\n");
		  fprintf(stderr,"Exiting...\n");
		  exit(1);
		}
		previous_sequence_length = sequence_length;

                for(i = 0; i < sequence_length/stride + 1; ++i){
                    for(j = 0; j < stride; ++j){
                        (GPU_BLASTP_database)[ base + (sequence + i*Group_size)*stride + j] = (seq_arg.seq->sequence)[AminoAcids];
                        AminoAcids++;
                        if( AminoAcids >= sequence_length )
                            goto done;
                    }
                }

              done:

                counter++;
                sequence++;

            } /* if ( seq_arg.oid < num_sequences_to ) */

            BlastSeqSrcReleaseSequence(seq_src, &seq_arg);

            if( seq_arg.oid == (break_limit-1) )
                break;

        } /* while loop */

        //Create the name of the gpu database
        if( 1 == num_volumes ) {
            int GPU_Database_name_bytes = strlen(Database_name) + strlen(".gpu") + 1;
            GPU_Database_name_temp = (char*) realloc(GPU_Database_name, GPU_Database_name_bytes*sizeof(char));
            if( NULL != GPU_Database_name_temp ){
                GPU_Database_name = GPU_Database_name_temp;
                memset(GPU_Database_name, 0, GPU_Database_name_bytes);
            } else {
                printf("Error (re)allocating memory for the GPU database creation\n");
                exit(1);
            }
            strcat(GPU_Database_name, Database_name);
        } else { //if there are more than one volumes then extra characters are needed for the volume number
            sprintf(digits,".%02d",volume);
            Int4 GPU_Database_name_bytes = strlen(Database_name) + strlen(digits) + strlen(".gpu") + 1;
            GPU_Database_name_temp = (char*) realloc(GPU_Database_name, GPU_Database_name_bytes*sizeof(char));
            if( NULL != GPU_Database_name_temp ){
                GPU_Database_name = GPU_Database_name_temp;
                memset(GPU_Database_name, 0, GPU_Database_name_bytes);
            } else {
                printf("Error (re)allocating memory for the GPU database creation\n");
                exit(1);
            }
            strcat(GPU_Database_name, Database_name);
            strcat(GPU_Database_name, digits);
        }

        strcat(GPU_Database_name, ".gpu");

        FILE* Database_file = fopen( GPU_Database_name,"wb");
        if ( NULL == Database_file )
            printf("GPU Database file cannot be opened for writing\n");

        if( GPU_VERBOSE )
            printf("GPU BLASTP database size = %ld bytes\n", GPU_BLASTP_database_bytes);

        fwrite(GPU_BLASTP_database, sizeof(char), GPU_BLASTP_database_bytes, Database_file);
        fclose(Database_file);

        Write_GPU_information_file(GPU_information_file,
                                   GPU_Database_name,
                                   GPU_BLASTP_database_bytes,
                                   volume_lengths[volume],
                                   num_sequences_to,
                                   store_limit,
                                   break_limit,
                                   Group_number,
                                   Sequence_length_vector);


        printf("Done with creating the GPU Database file (%s)\n", GPU_Database_name);

    }//for(volume = 0; volume < num_volumes; ++volume)

    fclose(GPU_information_file);
    free(GPU_Database_name);
    free(GPU_BLASTP_database);
    free(Sequence_length_vector);
    free(digits);

}

/** Read the database used by GPU_BLASTP
 * @param h_GPU_BLASTP_database        array that holds the GPU database [out]
 * @param h_GPU_BLASTP_database_bytes  size in bytes of the GPU_database [out]
*/

void ReadGPU_BLASTP_database(const char* h_GPU_BLASTP_database,
                             const char* GPU_volume_name,
                             const unsigned long h_GPU_BLASTP_database_bytes){

    FILE * pFile;
    long lSize;
    size_t result;

    pFile = fopen ( GPU_volume_name, "rb" );

    if (pFile==NULL) {fputs ("File error: GPU database file could not be oppened",stderr); exit (1);}

    lSize = h_GPU_BLASTP_database_bytes;

    result = fread (h_GPU_BLASTP_database, sizeof(char), lSize, pFile);

    fclose (pFile);

}


void GPU_BLASTP_Execute(const BlastSeqSrc* seq_src,
                        const LookupTableWrap* lookup_wrap,
                        const BlastCoreAuxStruct* aux_struct,
                        const BLAST_SequenceBlk* query,
                        const BlastQueryInfo* query_info,
                        const BlastInitialWordParameters* word_params,
                        const BlastGPUOptions* gpu_options,
                        const Int4 h_Hits_bytes,
                        const Int4 h_GPUBlastInitHitList_bytes,
                        const char* GPU_volume_name,
                        const unsigned long h_GPU_BLASTP_database_bytes,
                        const Int4* h_Sequence_length_vector,
                        const Int4 Group_number,
                        Int4** d_Hits,
                        GPUBlastInitHitList** d_GPUBlastInitHitList,
                        Int4** d_Database2Dpadded,
                        Int4** d_RepeatedSubstitutionMatrix,
                        Int4** d_RepeatedSequence_length_vector,
                        PV_ARRAY_TYPE** d_RepeatedPV,
                        AaLookupSmallboneCell** d_ThickBackbone,
                        Uint2** d_overflow ,
                        Int4** d_RepeatedDiag_array,
                        double *GPU_elapsed,
                        double *FillDatabase_elapsed)
{
    const char* h_GPU_BLASTP_database = (char*) calloc( h_GPU_BLASTP_database_bytes, sizeof(char) );

    struct timespec GPU_start, GPU_end, CPU_start, CPU_end;

    clock_gettime(CLOCK_MONOTONIC, &CPU_start);

    if(GPU_VERBOSE)
        fprintf(stderr,"threadId = 0: Reading the GPU database...\n");
    ReadGPU_BLASTP_database(h_GPU_BLASTP_database,
                            GPU_volume_name,
                            h_GPU_BLASTP_database_bytes);
    if(GPU_VERBOSE)
        fprintf(stderr,"threadId = 0: Done reading the GPU database\n");

    clock_gettime(CLOCK_MONOTONIC, &CPU_end);
    *FillDatabase_elapsed = (CPU_end.tv_sec - CPU_start.tv_sec) + (double)(CPU_end.tv_nsec - CPU_start.tv_nsec)/ BILLION;

    clock_gettime(CLOCK_MONOTONIC, &GPU_start);

    if( sufficient_memory )
        GPU_BLASTP(query, query_info, lookup_wrap, aux_struct->ewp, word_params,
                   gpu_options, h_Sequence_length_vector,
                   h_GPU_BLASTP_database_bytes, h_GPU_BLASTP_database,
                   d_Hits, h_Hits_bytes, d_GPUBlastInitHitList, h_GPUBlastInitHitList_bytes,
                   d_Database2Dpadded, d_RepeatedSubstitutionMatrix, d_RepeatedSequence_length_vector,
                   d_RepeatedPV, d_ThickBackbone, d_overflow, d_RepeatedDiag_array
                   );

    clock_gettime(CLOCK_MONOTONIC, &GPU_end);
    *GPU_elapsed = (double) (GPU_end.tv_sec - GPU_start.tv_sec) + (double)(GPU_end.tv_nsec - GPU_start.tv_nsec) / BILLION;

    free(h_GPU_BLASTP_database);

}

/** Perform the gapped alignment of the hits returned by the GPU
 * @param  program_number Type of BLAST program [in]
 * @param query The query sequence [in]
 * @param query_info Additional query information [in]
 * @param seq_src Structure containing BLAST database [in]
 * @param gap_align Structure containing scoring block and memory allocated
 *                  for gapped alignment. [in]
 * @param score_params Hit scoring parameters [in]
 * @param lookup_wrap The lookup table, constructed earlier [in]
 * @param word_options Options for processing initial word hits [in]
 * @param ext_params Parameters for the gapped extension [in]
 * @param hit_params Parameters for saving the HSPs [in]
 * @param eff_len_params Parameters for setting effective lengths [in]
 * @param psi_options Options specific to PSI-BLAST [in]
 * @param db_options Options for handling BLAST database [in]
 * @param gpu_options Options for handling the GPU execution [in]
 * @param hsp_stream Placeholder for saving HSP lists [in]
 * @param diagnostics Return statistics containing numbers of hits on
 *                    different stages of the search. Statistics saved only
 *                    for the allocated parts of the structure. [in] [out]
 * @param interrupt_search function callback to allow interruption of BLAST
 * search [in, optional]
 * @param progress_info contains information about the progress of the current
 * BLAST search [in|out]
 * @param h_Hits vector holding the hits, successful hits, and extensions returned by the GPU [in]
 * @param h_GPUBlastInitHitList vector holding the info of the successful hits returnedby the GPU [in]
 */


Int4
BLASTP_GetGPUHits(EBlastProgramType program_number,
                  BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
                  const BlastSeqSrc* seq_src, BlastGapAlignStruct* gap_align,
                  BlastScoringParameters* score_params,
                  LookupTableWrap* lookup_wrap,
                  const BlastInitialWordOptions* word_options,
                  BlastExtensionParameters* ext_params,
                  BlastHitSavingParameters* hit_params,
                  BlastEffectiveLengthsParameters* eff_len_params,
                  const PSIBlastOptions* psi_options,
                  const BlastDatabaseOptions* db_options,
                  BlastGPUOptions* gpu_options,
                  BlastHSPStream* hsp_stream, BlastDiagnostics* diagnostics,
                  TInterruptFnPtr interrupt_search, SBlastProgress* progress_info,
                  Int4* h_Hits, GPUBlastInitHitList* h_GPUBlastInitHitList,
                  BlastHSPList** hsp_list, BlastCoreAuxStruct* aux_struct,
                  BlastInitialWordParameters* word_params,
                  Int4 gpu_limit,
                  Int4 break_limit,
                  BlastSeqSrcIterator* itr,
                  BlastSeqSrcGetSeqArg seq_arg,
                  Int4 volume_base,
                  Int4 threadId){

    Int2 status = 0;
    Int8 db_length = 0;
    const BlastScoringOptions* score_options = score_params->options;
    Boolean gapped_calculation = score_options->gapped_calculation;
    BlastScoreBlk* sbp = gap_align->sbp;

    const Boolean kNucleotide = (program_number == eBlastTypeBlastn ||
                                program_number == eBlastTypePhiBlastn);

   /* remember the current search state */
   if (progress_info)
       progress_info->stage = ePrelimSearch;

   /* For RPS BLAST, there is no loop over subject sequences, so the preliminary
      search engine is done in a separate function. */

   /* Update the parameters for linking HSPs, if necessary. */



   /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
      sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
   seq_arg.encoding = eBlastEncodingProtein;

   db_length = BlastSeqSrcGetTotLen(seq_src);

   Int4 print_point = BlastSeqSrcGetNumSeqs(seq_src) / 50;

/* iterate over all subject sequences */
   while ( ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) ) {
       Int4 stat_length;
       if (seq_arg.oid == BLAST_SEQSRC_ERROR)
           break;
       if (BlastSeqSrcGetSequence(seq_src, (void*) &seq_arg) < 0)
           continue;
       if (db_length == 0) {
           /* This is not a database search, hence need to recalculate and save
              the effective search spaces and length adjustments for all
              queries based on the length of the current single subject
              sequence. */
           if ((status = BLAST_OneSubjectUpdateParameters(program_number,
                                                          seq_arg.seq->length, score_options, query_info,
                                                          sbp, hit_params, word_params,
                                                          eff_len_params)) != 0)
               return status;
       }

       stat_length = seq_arg.seq->length;

       /* Calculate cutoff scores for linking HSPs. Do this only for
          ungapped protein searches and ungapped translated searches. */
       if (hit_params->link_hsp_params && !kNucleotide && !gapped_calculation) {
           CalculateLinkHSPCutoffs(program_number, query_info, sbp,
                                   hit_params->link_hsp_params, word_params, db_length,
                                   seq_arg.seq->length);
       }

       if (Blast_SubjectIsTranslated(program_number)) {
           /* If the subject is translated and the BlastSeqSrc implementation
            * doesn't provide a genetic code string, use the default genetic
            * code for all subjects (as in the C toolkit) */
           if (seq_arg.seq->gen_code_string == NULL) {
               seq_arg.seq->gen_code_string =
                   GenCodeSingletonFind(db_options->genetic_code);
           }
           ASSERT(seq_arg.seq->gen_code_string);
           stat_length /= CODON_LENGTH;                    
       }
       
       /* If the subject is translated and the BlastSeqSrc implementation
        * doesn't provide a genetic code string, use the default genetic
        * code for all subjects (as in the C toolkit) */
       if( seq_arg.oid < gpu_limit ){

           if( GPU_VERBOSE )
               if( !(seq_arg.oid % print_point ) )
                   fprintf(stderr,"threadId = %d, seq_arg.oid = %6d, ranges:[%6d, %6d)\n",
                           threadId, seq_arg.oid, itr->oid_range[0], itr->oid_range[1]);

           status =
               s_BlastSearchEngineCore(program_number, query, query_info,
                                       seq_arg.seq, lookup_wrap, gap_align, score_params, word_params,
                                       ext_params, hit_params, db_options, diagnostics, aux_struct,
                                       &(*hsp_list), interrupt_search, progress_info,
                                       gpu_options, h_Hits, h_GPUBlastInitHitList, volume_base);
       }

       if (status) {
           break;
       }

       if (*hsp_list && (*hsp_list)->hspcnt > 0) {
           int query_index=0; /* Used to loop over queries below. */
           if (!gapped_calculation) {
               /* The following must be performed for any ungapped
                  search with a nucleotide database. */
               
               status =
                   Blast_HSPListReevaluateUngapped(
                       program_number, *hsp_list, query,
                       seq_arg.seq, word_params, hit_params,
                       query_info, sbp, score_params, seq_src,
                       seq_arg.seq->gen_code_string);
               
               if (status) {
                   BlastSeqSrcReleaseSequence(seq_src, (void*) &seq_arg);
                   return status;
               }
               /* Relink HSPs if sum statistics is used, because scores might
                * have changed after reevaluation with ambiguities, and there
                * will be no traceback stage where relinking is done normally.
                * If sum statistics are not used, just recalculate e-values.
                */
               if (hit_params->link_hsp_params) {
                   status =
                       BLAST_LinkHsps(program_number, *hsp_list, query_info,
                                      seq_arg.seq->length, sbp,
                                      hit_params->link_hsp_params,
                                      gapped_calculation);
               } else {
//                   Blast_HSPListGetEvalues(query_info, stat_length,
//                                           hsp_list, gapped_calculation, FALSE,
//                                           sbp, 0, 1.0);
                   Blast_HSPListGetEvalues(program_number, query_info,
                                           seq_arg.seq->length, *hsp_list,
                                           gapped_calculation, FALSE,
                                           sbp, 0, 1.0);
               }
               /* Use score threshold rather than evalue if 
                * matrix_only_scoring is used.  -RMH- 
                */
               if ( sbp->matrix_only_scoring )
               {
                   status = Blast_HSPListReapByRawScore(*hsp_list,
                                                        hit_params->options);
               }else {
                   status = Blast_HSPListReapByEvalue(*hsp_list,
                                                      hit_params->options);
               }
               
               /* Calculate and fill the bit scores, since there will be no
                  traceback stage where this can be done. */
               Blast_HSPListGetBitScores(*hsp_list, gapped_calculation, sbp);
           } /* if (!gapped_calculation) */
           
           /* Save the results. */
           status = BlastHSPStreamWrite(hsp_stream, hsp_list);
           if (status != 0)
               break;

           if (hit_params->low_score)
           {
               for (query_index=0; query_index<hsp_stream->results->num_queries; query_index++)
                   if (hsp_stream->results->hitlist_array[query_index] && hsp_stream->results->hitlist_array[query_index]->heapified)
                       hit_params->low_score[query_index] = 
                           MAX(hit_params->low_score[query_index], 
                               hit_params->options->low_score_perc*(hsp_stream->results->hitlist_array[query_index]->low_score));
           }
       } /* if (hsp_list && hsp_list->hspcnt > 0) */
       
       
       BlastSeqSrcReleaseSequence(seq_src, (void*) &seq_arg);
       
       /* check for interrupt */
       if (interrupt_search && (*interrupt_search)(progress_info) == TRUE) {
           status = BLASTERR_INTERRUPTED;
           break;
       }
       
       if( seq_arg.oid >= (break_limit-1) )
           break;
       
   }//while loop
   
   /* Fill the cutoff values in the diagnostics structure */
   if (diagnostics && diagnostics->cutoffs) {
       s_FillReturnCutoffsInfo(diagnostics->cutoffs, score_params, word_params, ext_params, hit_params);
   }

   return status;
   
}



Int4
GPU_BLAST_PreliminarySearchEngine(EBlastProgramType program_number,
                                  BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
                                  const BlastSeqSrc* seq_src, BlastGapAlignStruct* gap_align,
                                  BlastScoringParameters* score_params,
                                  LookupTableWrap* lookup_wrap,
                                  const BlastInitialWordOptions* word_options,
                                  BlastExtensionParameters* ext_params,
                                  BlastHitSavingParameters* hit_params,
                                  BlastEffectiveLengthsParameters* eff_len_params,
                                  const PSIBlastOptions* psi_options,
                                  const BlastDatabaseOptions* db_options,
                                  const BlastGPUOptions* gpu_options,
                                  BlastHSPStream* hsp_stream, BlastDiagnostics* diagnostics,
                                  TInterruptFnPtr interrupt_search, SBlastProgress* progress_info)
{
    BlastCoreAuxStruct* aux_struct = NULL;
    BlastHSPList* hsp_list = NULL;
    BlastSeqSrcGetSeqArg seq_arg;
    BlastSeqSrcGetSeqArg gpu_seq_arg;
    Int2 status = 0;
    Int8 db_length = 0;
    const BlastScoringOptions* score_options = score_params->options;
    const BlastHitSavingOptions* hit_options = hit_params->options;
    const BlastExtensionOptions* ext_options = ext_params->options;
    BlastInitialWordParameters* word_params = NULL;
    Boolean gapped_calculation = score_options->gapped_calculation;
    BlastScoreBlk* sbp = gap_align->sbp;
    BlastSeqSrcIterator* itr;
    BlastSeqSrcIterator* gpu_itr;

    itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/2000,1));
    gpu_itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/2000,1));

    memset((void*) &seq_arg, 0, sizeof(seq_arg));
    memset((void*) &gpu_seq_arg, 0, sizeof(gpu_seq_arg));


    const Boolean kNucleotide = (program_number == eBlastTypeBlastn ||
                                 program_number == eBlastTypePhiBlastn);

    BlastInitialWordParametersNew(program_number, word_options,
                                  hit_params, lookup_wrap, sbp, query_info,
                                  BlastSeqSrcGetAvgSeqLen(seq_src), &word_params);

    if ((status =
         s_BlastSetUpAuxStructures(seq_src, lookup_wrap, word_params,
                                   ext_options, hit_options, query, &aux_struct)) != 0)
        return status;

    pthread_mutex_lock(&thread_counter_mutex);
    Int4 threadId = thread_counter++;
    if( 0 == threadId ){
        int init_return = pthread_barrier_init(&barr, NULL, gpu_options->cpu_threads);
    }
    pthread_mutex_unlock(&thread_counter_mutex);

    if ( (2 == gpu_options->method) ){
        if( 0 == threadId ) {//only one threads creates the database
            CreateGPU_BLASTP_database(seq_src, gpu_options);
        }
	exit(0);
    }

    struct timespec CPU_start, CPU_end, Total_start, Total_end;

    double GPU_elapsed = 0, CPU_elapsed = 0, Total_elapsed = 0, FillDatabase_elapsed = 0;

    Int4 num_blocksx, num_threadsx, num_volumes = 0, percentage = 0;
    FILE* GPU_information_file = NULL;

    Read_GPU_information_file_preliminary(&GPU_information_file, seq_src,
                                          &num_blocksx, &num_threadsx,
                                          &percentage, &num_volumes);

    //By using gpu_options_local there is no need for a mutex and synchronization for writing "gpu_options_local->use_gpu = FALSE" later on,
    BlastGPUOptions* gpu_options_local = (BlastGPUOptions *) calloc(1, sizeof(BlastGPUOptions));
    gpu_options_local->use_gpu = gpu_options->use_gpu;
    gpu_options_local->num_blocksx = num_blocksx;
    gpu_options_local->num_threadsx = num_threadsx;
    gpu_options_local->num_sequences_to = percentage;
    gpu_options_local->method = gpu_options->method;
    char* GPU_volume_name = NULL;
    Int4* Sequence_length_vector = NULL;
    Int4* h_Hits_temp = NULL;
    GPUBlastInitHitList* h_GPUBlastInitHitList_temp = NULL;

    //Initialize global variables
    if( 0 == threadId ) {
      h_Hits = NULL;
      h_GPUBlastInitHitList = NULL;
      sufficient_memory = FALSE;
    }

    Int4 volume = 0;
    Int4 volume_length = 0;
    Int4 gpu_NextChunkOID = 0;
    Int4 break_limit = 0;
    for(volume = 0; volume < num_volumes; ++volume)
    {
        clock_gettime(CLOCK_MONOTONIC, &Total_start);
        unsigned long h_GPU_BLASTP_database_bytes = 0;
        Int4 num_sequences_to = 0;
        Int4 gpu_limit = 0, gpu_limit_temp = 0;
	Int4 break_limit_temp = 0;
        Int4 Group_number = 0;
        gpu_NextChunkOID += volume_length;

	Read_GPU_information_file(GPU_information_file,
                                  &GPU_volume_name,
                                  &h_GPU_BLASTP_database_bytes, &volume_length,
				  &num_sequences_to, &gpu_limit_temp,
				  //&num_sequences_to, &gpu_limit,
				  &break_limit_temp, &Group_number,
				  //&break_limit, &Group_number,
                                  &Sequence_length_vector,
				  gpu_options->cpu_threads);

	gpu_options_local->num_sequences_to = num_sequences_to;
	gpu_limit = break_limit + num_sequences_to;
	break_limit = break_limit_temp;

        /*num_sequences_to is multiplied by 3 because for each sequence we store hits, successful_extensions, hits_extended */
        Int4 h_Hits_bytes = num_sequences_to * 3 * sizeof(Int4);
	if( 0 == threadId ) {
	  h_Hits_temp = (Int4*) realloc(h_Hits_temp, h_Hits_bytes);
	  if( NULL != h_Hits_temp )
	    h_Hits = h_Hits_temp;
	  else {
	    puts("Error (re)allocating memory for the GPU\n");
	    exit(1);
	  }
	}

        Int4 h_GPUBlastInitHitList_bytes = num_sequences_to * MAX_SUCCESSFUL_HITS_PER_SEQUENCE * sizeof(GPUBlastInitHitList);
	if( 0 == threadId ) {
	  h_GPUBlastInitHitList_temp = (GPUBlastInitHitList*) realloc(h_GPUBlastInitHitList, h_GPUBlastInitHitList_bytes);
	  if( NULL != h_GPUBlastInitHitList_temp )
            h_GPUBlastInitHitList = h_GPUBlastInitHitList_temp;
	  else {
            printf("Error (re)allocating memory for the GPU\n");
            exit(1);
	  }
	}

        //The GPU arrays are declared here because we want to be able to deallocate them in this function. If we deallocated them in
        //GPU_BLASTP() then that would act as an thread synchronizer, and wouldn't allow asynchronous execution
        Int4 *d_Hits = NULL;
        GPUBlastInitHitList* d_GPUBlastInitHitList = NULL;
        Int4* d_Database2Dpadded = NULL;
        Int4* d_RepeatedSubstitutionMatrix = NULL;
        Int4* d_RepeatedSequence_length_vector = NULL;
        PV_ARRAY_TYPE* d_RepeatedPV = NULL;
        AaLookupSmallboneCell* d_ThickBackbone = NULL;
        Uint2* d_overflow = NULL;
        Int4* d_RepeatedDiag_array = NULL;

        sufficient_memory = GPU_BLASTP_check_memory(lookup_wrap, aux_struct->ewp,
                                                    gpu_options_local,
                                                    h_GPU_BLASTP_database_bytes,
                                                    h_Hits_bytes, h_GPUBlastInitHitList_bytes,
                                                    Group_number, query_info->num_queries,
                                                    query->length
                                                    );

	if( (0 == threadId) && !sufficient_memory ){
	  fprintf(stderr,"WARNING: Not enough GPU global memory to process volume No. %02d of the database. Continuing without the CPU...\n",volume);
	  fprintf(stderr,"         Consider splitting the input database in volumes with smaller size by using the option \"-max_file_size <String>\" (e.g. -max_file_sz 500MB).\n");
	  fprintf(stderr,"         Consider using fewer GPU blocks and then fewer GPU thread when formating the database with (e.g. \"-gpu_blocks 256\" and/or \"-gpu_threads 32\") to reduce the GPU global memory requiremenets.\n");
	}

	if( (0 == threadId) && sufficient_memory && (gpu_options->cpu_threads > 1) )
 	  BlastSeqSrcSetMmapMemory(seq_src, FALSE);
	else if ( (0 == threadId) && sufficient_memory  )
	  BlastSeqSrcSetMmapMemory(seq_src, TRUE);

        //wait until GPU_BLASTP_check_memory() is finished in order to have global variable
        //"sufficient_memory" updated for all threads. And, until MmapMemory is set
        pthread_barrier_wait(&barr);

        if( sufficient_memory && (0 == threadId) ) {
	  GPU_BLASTP_Execute(seq_src, lookup_wrap, aux_struct, query, query_info, word_params,
			     gpu_options_local,  h_Hits_bytes, h_GPUBlastInitHitList_bytes,
			     GPU_volume_name, h_GPU_BLASTP_database_bytes,
			     Sequence_length_vector, Group_number,
			     &d_Hits, &d_GPUBlastInitHitList, &d_Database2Dpadded,
			     &d_RepeatedSubstitutionMatrix, &d_RepeatedSequence_length_vector,
			     &d_RepeatedPV, &d_ThickBackbone, &d_overflow , &d_RepeatedDiag_array,
			     &GPU_elapsed, &FillDatabase_elapsed);
	}

        /* remember the current search state */
        if (progress_info)
	  progress_info->stage = ePrelimSearch;

        /* For RPS BLAST, there is no loop over subject sequences, so the preliminary
           search engine is done in a separate function. */
        if (Blast_ProgramIsRpsBlast(program_number)) {
	  status =
	    s_RPSPreliminarySearchEngine(program_number, query, query_info,
					 seq_src, score_params, lookup_wrap, aux_struct, word_params,
					 ext_params, gap_align, hit_params, hsp_stream, diagnostics,
					 interrupt_search, progress_info,
					 gpu_options_local, h_Hits, h_GPUBlastInitHitList );
	  word_params = BlastInitialWordParametersFree(word_params);
	  s_BlastCoreAuxStructFree(aux_struct);
	  return status;
        }

        /* Update the parameters for linking HSPs, if necessary. */
        BlastLinkHSPParametersUpdate(word_params, hit_params, gapped_calculation);

        /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
           sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
        seq_arg.encoding = eBlastEncodingProtein;

        db_length = BlastSeqSrcGetTotLen(seq_src);

        //gpu_options_local->use_gpu is set to FALSE in order the CPU to carry out the rest of alignments
        gpu_options_local->use_gpu = FALSE;
        //if there is not enough GPU memory use only the CPU and start alignment from the first database sequence
        num_sequences_to = sufficient_memory ? num_sequences_to : 0;
        gpu_limit = sufficient_memory ? gpu_limit : (break_limit - volume_length);

	Int4 print_point = 0;
        if(GPU_VERBOSE){
            fprintf(stderr,"threadID = %d: CPU execution...\n", threadId);
	    print_point = BlastSeqSrcGetNumSeqs(seq_src) / 50;
	}

        /* iterate over all subject sequences */
        clock_gettime(CLOCK_MONOTONIC, &CPU_start);
        while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
            Int4 stat_length;
            if (seq_arg.oid == BLAST_SEQSRC_ERROR)
                break;
            
            if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0)
                continue;
            
            if( seq_arg.oid >= gpu_limit  ){

                if(GPU_VERBOSE)
                    if( !(seq_arg.oid % print_point ) )
                        fprintf(stderr,"threadId = %d, seq_arg.oid = %6d, ranges:[%6d, %6d)\n",
                                threadId, seq_arg.oid, itr->oid_range[0], itr->oid_range[1]);

                if (db_length == 0) {
                    /* This is not a database search, hence need to recalculate and save
                       the effective search spaces and length adjustments for all
                       queries based on the length of the current single subject
                       sequence. */
                    if ((status = BLAST_OneSubjectUpdateParameters(program_number,
                                                                   seq_arg.seq->length, score_options, query_info,
                                                                   sbp, hit_params, word_params,
                                                                   eff_len_params)) != 0)
                        return status;
                }

                stat_length = seq_arg.seq->length;
                
                /* Calculate cutoff scores for linking HSPs. Do this only for
                   ungapped protein searches and ungapped translated
                   searches. */
                if (hit_params->link_hsp_params && !kNucleotide && !gapped_calculation) {
                    CalculateLinkHSPCutoffs(program_number, query_info, sbp,
                                            hit_params->link_hsp_params, word_params, db_length,
                                            seq_arg.seq->length);
                }


                if (Blast_SubjectIsTranslated(program_number)) {
                    /* If the subject is translated and the BlastSeqSrc implementation
                     * doesn't provide a genetic code string, use the default genetic
                     * code for all subjects (as in the C toolkit) */
                    if (seq_arg.seq->gen_code_string == NULL) {
                        seq_arg.seq->gen_code_string =
                            GenCodeSingletonFind(db_options->genetic_code);
                    }
                    ASSERT(seq_arg.seq->gen_code_string);
                    stat_length /= CODON_LENGTH;                    
                }

                status =
                    s_BlastSearchEngineCore(program_number, query, query_info,
                                            seq_arg.seq, lookup_wrap, gap_align, score_params, word_params,
                                            ext_params, hit_params, db_options, diagnostics, aux_struct,
                                            &hsp_list, interrupt_search, progress_info,
                                            gpu_options_local, h_Hits, h_GPUBlastInitHitList, 0);

                if (status) {
                    break;
                }

                if (hsp_list && hsp_list->hspcnt > 0) {
                    int query_index=0; /* Used to loop over queries below. */
                    if (!gapped_calculation) {
                        /* The following must be performed for any ungapped
                           search with a nucleotide database. */
                        status =
                            Blast_HSPListReevaluateUngapped(
                                program_number, hsp_list, query,
                                seq_arg.seq, word_params, hit_params,
                                query_info, sbp, score_params, seq_src,
                                seq_arg.seq->gen_code_string);
                        
                        if (status) {
                            BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
                            return status;
                        }
                        /* Relink HSPs if sum statistics is used, because scores might
                         * have changed after reevaluation with ambiguities, and there
                         * will be no traceback stage where relinking is done normally.
                         * If sum statistics are not used, just recalculate e-values.
                         */
                        if (hit_params->link_hsp_params) {
                            status =
                                BLAST_LinkHsps(program_number, hsp_list, query_info,
                                               seq_arg.seq->length, sbp,
                                               hit_params->link_hsp_params,
                                               gapped_calculation);
                        } else {
//                            Blast_HSPListGetEvalues(query_info, stat_length,
//                                                    hsp_list, gapped_calculation, FALSE,
//                                                    sbp, 0, 1.0);
                            Blast_HSPListGetEvalues(program_number,query_info,
                                                    seq_arg.seq->length, hsp_list,
                                                    gapped_calculation,
                                                    0, sbp, 0, 1.0);
                        }
                        /* Use score threshold rather than evalue if 
                         * matrix_only_scoring is used.  -RMH- 
                         */
                        if ( sbp->matrix_only_scoring )
                        {
                            status = Blast_HSPListReapByRawScore(hsp_list,
                                                                 hit_params->options);
                        }else {
                            status = Blast_HSPListReapByEvalue(hsp_list,
                                                               hit_params->options);
                        }
                        
                        /* Calculate and fill the bit scores, since there will be no
                           traceback stage where this can be done. */
                        Blast_HSPListGetBitScores(hsp_list, gapped_calculation, sbp);
                    }
                    
                    /* Save the results. */
                    status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
                    if (status != 0)
                        break;

                    if (hit_params->low_score)
                    {
                        for (query_index=0; query_index<hsp_stream->results->num_queries; query_index++)
                            if (hsp_stream->results->hitlist_array[query_index] && hsp_stream->results->hitlist_array[query_index]->heapified)
                                hit_params->low_score[query_index] = 
                                    MAX(hit_params->low_score[query_index], 
                                        hit_params->options->low_score_perc*(hsp_stream->results->hitlist_array[query_index]->low_score));
                    }

                }//if (hsp_list && hsp_list->hspcnt > 0) 
            } //if( seq_arg.oid >= gpu_limit  )

            BlastSeqSrcReleaseSequence(seq_src, &seq_arg);

            /* check for interrupt */
            if (interrupt_search && (*interrupt_search)(progress_info) == TRUE) {
                status = BLASTERR_INTERRUPTED;
                break;
            }

            if( seq_arg.oid >= (break_limit-1) )
                break;
        }// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) 

        clock_gettime(CLOCK_MONOTONIC, &CPU_end);
        CPU_elapsed = (CPU_end.tv_sec - CPU_start.tv_sec) + (double)(CPU_end.tv_nsec - CPU_start.tv_nsec) / BILLION;

        if(GPU_VERBOSE)
	  fprintf(stderr,"threadID = %d: Done with CPU processing (%.2f sec)\n", threadId, CPU_elapsed);

        clock_gettime(CLOCK_MONOTONIC, &CPU_start);

        //writing h_Hits and h_GPUBlastInitHitList doesn't require a mutex because because it is accessed by only one thread
        if( sufficient_memory && (0 == threadId) )
            GPU_BLASTP_get_data(h_Hits, d_Hits, h_Hits_bytes,
                                h_GPUBlastInitHitList, d_GPUBlastInitHitList, h_GPUBlastInitHitList_bytes,
                                num_sequences_to);

        //NOTE: GPU_BLASTP_get_data() acts as a synchronizer for the GPU since it calls cudaMemcpy(); i.e., GPU_BLASTP_get_data_elapsed possibly
        //includes some waiting time for the GPU to finish processing
        clock_gettime(CLOCK_MONOTONIC, &CPU_end);
        double GPU_BLASTP_get_data_elapsed = (CPU_end.tv_sec - CPU_start.tv_sec) + (double)(CPU_end.tv_nsec - CPU_start.tv_nsec) / BILLION;


        clock_gettime(CLOCK_MONOTONIC, &CPU_start);
        //wait for both threads before calling BlastSeqSrcResetChunkIterator() since one thread might reset seq_src while the
        //other one still working on the database
        if( sufficient_memory )
            pthread_barrier_wait(&barr);

        clock_gettime(CLOCK_MONOTONIC, &CPU_end );
        double waiting_time1 = ( CPU_end.tv_sec - CPU_start.tv_sec ) + (double)( CPU_end.tv_nsec - CPU_start.tv_nsec ) / BILLION;

        clock_gettime(CLOCK_MONOTONIC, &CPU_start);

        Int4 cpu_NextChunkOID = itr->oid_range[1];
        if( sufficient_memory && (0 == threadId) )
            BlastSeqSrcSetChunkIterator(seq_src, gpu_NextChunkOID);

        clock_gettime(CLOCK_MONOTONIC, &CPU_start);
        //wait for both threads after BlastSeqSrcResetChunkIterator() since one thread might start working on the database
        //in BLASTP_GetGPUHits() and it can find seq_src pointing at the end of the database, and thus won't work on it.
        pthread_barrier_wait(&barr);

        clock_gettime(CLOCK_MONOTONIC, &CPU_end );
        double waiting_time2 = ( CPU_end.tv_sec - CPU_start.tv_sec ) + (double)( CPU_end.tv_nsec - CPU_start.tv_nsec ) / BILLION;

        double BLAST_GetGPUHits_elapsed = 0;
        if( sufficient_memory )
        {
            clock_gettime(CLOCK_MONOTONIC, &CPU_start );

            if(GPU_VERBOSE)
                fprintf(stderr,"threadID = %d: CPU processing GPU output...\n",threadId);

            gpu_options_local->use_gpu = TRUE;
            gpu_itr->current_pos = UINT4_MAX;//make this UINT4_MAX to get a new chunk

	    //When processing the GPU outputs use mmap since all sequences are processed by the CPU
	    if( 0 == threadId )
	      BlastSeqSrcSetMmapMemory(seq_src, FALSE );

            status = BLASTP_GetGPUHits(program_number, query, query_info, seq_src, gap_align,
                                       score_params, lookup_wrap, word_options, ext_params,
                                       hit_params, eff_len_params, psi_options, db_options,
                                       gpu_options_local, hsp_stream, diagnostics, interrupt_search,
                                       progress_info, h_Hits, h_GPUBlastInitHitList,
                                       &hsp_list, aux_struct, word_params, gpu_limit, break_limit,
                                       gpu_itr, gpu_seq_arg, gpu_NextChunkOID,
                                       threadId);

            clock_gettime(CLOCK_MONOTONIC, &CPU_end );
            BLAST_GetGPUHits_elapsed = ( CPU_end.tv_sec - CPU_start.tv_sec ) + (double)( CPU_end.tv_nsec - CPU_start.tv_nsec ) / BILLION;

            if(GPU_VERBOSE)
                fprintf(stderr,"threadId = %d: Done with CPU processing GPU output (%.2f sec)\n",threadId, BLAST_GetGPUHits_elapsed);

        }

        clock_gettime(CLOCK_MONOTONIC, &CPU_start);

        //wait until both threads are done with BLASTP_GetGPUHits() before deallocating h_Hits and h_GPUBlastInitHitList since
        //if one threadId=0 deallocates these arrays while other threads are still working on them a segmentation fault will occur
        if( sufficient_memory )
            pthread_barrier_wait(&barr);

        clock_gettime(CLOCK_MONOTONIC, &CPU_end );
        double waiting_time3 = ( CPU_end.tv_sec - CPU_start.tv_sec ) + (double)( CPU_end.tv_nsec - CPU_start.tv_nsec ) / BILLION;

        clock_gettime(CLOCK_MONOTONIC, &CPU_start);

        if( sufficient_memory && (0 == threadId) ){
            BlastSeqSrcSetChunkIterator(seq_src, cpu_NextChunkOID);
            GPU_BLASTP_free_memory( &d_Hits, &d_GPUBlastInitHitList,
                                    &d_Database2Dpadded, &d_RepeatedSubstitutionMatrix, &d_RepeatedSequence_length_vector,
                                    &d_RepeatedPV, &d_ThickBackbone, &d_overflow, &d_RepeatedDiag_array);

        }

        clock_gettime(CLOCK_MONOTONIC, &CPU_end );
        double freeing_memory = ( CPU_end.tv_sec - CPU_start.tv_sec ) + (double)( CPU_end.tv_nsec - CPU_start.tv_nsec ) / BILLION;

        clock_gettime(CLOCK_MONOTONIC, &Total_end);
        Total_elapsed = (Total_end.tv_sec - Total_start.tv_sec) + (double)(Total_end.tv_nsec - Total_start.tv_nsec)/ BILLION;


        if(GPU_VERBOSE) {
            pthread_mutex_lock(&thread_print_mutex);
            fprintf(stderr,"\nthreadId = %d: Fill GPU database                   = %.2f sec",threadId, FillDatabase_elapsed);
            fprintf(stderr,"\nthreadId = %d: GPU processing                      = %.2f sec",threadId,GPU_elapsed);
            fprintf(stderr,"\nthreadId = %d: CPU processing                      = %.2f sec",threadId,CPU_elapsed);
            fprintf(stderr,"\nthreadId = %d: Get GPU data                        = %.2f sec",threadId,GPU_BLASTP_get_data_elapsed);
            fprintf(stderr,"\nthreadId = %d: Waiting before iterator reset       = %.2f sec",threadId,waiting_time1);
            fprintf(stderr,"\nthreadId = %d: Waiting after iterator reset        = %.2f sec",threadId,waiting_time2);
            fprintf(stderr,"\nthreadId = %d: GPU output processing               = %.2f sec",threadId,BLAST_GetGPUHits_elapsed);
            fprintf(stderr,"\nthreadId = %d: Waiting after GPU output processing = %.2f sec",threadId,waiting_time3);
            fprintf(stderr,"\nthreadId = %d: Freeing GPU memory                  = %.2f sec",threadId,freeing_memory);
            fprintf(stderr,"\nthreadId = %d: Total                               = %.2f sec\n\n",threadId,Total_elapsed);
            pthread_mutex_unlock(&thread_print_mutex);
        }
    } //for(volume = 0; volume < num_volumes; ++volume)

    if( sufficient_memory && (0 == threadId) ){
        free(h_Hits);
        free(h_GPUBlastInitHitList);
    }

    pthread_mutex_lock(&thread_counter_mutex);
    thread_counter--;
    pthread_mutex_unlock(&thread_counter_mutex);

    free(Sequence_length_vector);
    free(GPU_volume_name);
    free(gpu_options_local);
    fclose(GPU_information_file);


    hsp_list = Blast_HSPListFree(hsp_list);  /* in case we were interrupted */
    BlastSequenceBlkFree(seq_arg.seq);
    itr = BlastSeqSrcIteratorFree(itr);
    gpu_itr = BlastSeqSrcIteratorFree(gpu_itr);

    /* Fill the cutoff values in the diagnostics structure */
    if (diagnostics && diagnostics->cutoffs) {
        s_FillReturnCutoffsInfo(diagnostics->cutoffs, score_params, word_params,
                                ext_params, hit_params);
    }

    word_params = BlastInitialWordParametersFree(word_params);
    s_BlastCoreAuxStructFree(aux_struct);

    return status;
}
