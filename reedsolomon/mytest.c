#include <string.h>
#include <time.h>

#include "Headers/op.h"
#include "Headers/array.h"
#include "Headers/encode.h"
#include "Headers/decode.h"

int main()
{
	int BUFFER_SIZE = 255-8; //100; //max number of bytes to encode, in ASCII, deleted if over BUFFER_SIZE
	int ENCODE_SYMBOLS = 8; //encode bytes, BUFFER_SIZE + ENCODE_SYMBOLS <= 255
	int ERRORS = 4; //correctable errors, in bytes, <= ENCODE_SYMBOLS/2
	double start_time_encode, end_time_encode, start_time_decode, end_time_decode;
	
	time_t t;
	srand((unsigned) time(&t));

	struct gf_tables *gf_table = malloc(sizeof(struct gf_tables));
	gf_table->gf_exp = malloc(sizeof(struct Array));
	gf_table->gf_log = malloc(sizeof(struct Array));
	initArray(gf_table->gf_exp, 512);
	initArray(gf_table->gf_log, 256);
	gf_table = init_tables();
	
	printf("##### Reed-Solomon Error Correction #####\n");
	// printf("Enter the string you want to test and press enter when you're done:\n");
	
	struct Array *msg_in = malloc(sizeof(struct Array));
	initArray(msg_in, 50); //BUFFER_SIZE/2
	char my_msg[BUFFER_SIZE];
	
	//fgets( my_msg, sizeof(my_msg), stdin );
	for(int i = 0; i < BUFFER_SIZE; i++){
		my_msg[i] = 48 + rand() % 10;
	}

	for (size_t i = 0; i < strlen(my_msg); i++) {
		msg_in->array[i] = (int)my_msg[i];
		insertArray(msg_in); //msg_in->used++, and msg_in->size *= 2 if size not enough
	}
	printf("]\n");

	printf("Msg in: [");
	for (size_t i = 0; i < msg_in->used; i++) {
		printf("%u,",msg_in->array[i]);
	}
	printf("]\n");

	struct Array *msg = malloc(sizeof(struct Array)); //msg_in + encode
	initArray(msg, 170);

	struct Array *err_loc = malloc(sizeof(struct Array));

	struct Array *synd = malloc(sizeof(struct Array));

	struct Array *pos = malloc(sizeof(struct Array));

	struct Array *rev_pos = malloc(sizeof(struct Array));

	start_time_encode = clock();
	msg = rs_encode_msg(msg_in, ENCODE_SYMBOLS, gf_table); 
	end_time_encode = clock();
	printf("encode time=%f\n", (double)(end_time_encode-start_time_encode)/CLOCKS_PER_SEC);  
	
	printf("Msg Encoded: [");
	for (size_t i = 0; i < msg->used; i++) {
		printf("%u,",msg->array[i]);
	}
	printf("]\n");
	
	
	//Tempering msg
	// msg->array[0] = '0';
	// msg->array[3] = '0';
	// msg->array[10] = '0';
	for (int i = 0; i < ERRORS*3; i=i+3){ //random
		msg->array[i] = '0';
	}
	
	printf("Msg Tempered: [");
	for (size_t i = 0; i < strlen(my_msg); i++) {
		printf("%u,",msg->array[i]);
	}
	printf("]\n");
	
	printf("Msg Tempered: ");
	for (size_t i = 0; i < strlen(my_msg); i++) {
		printf("%c",msg->array[i]);
	}
	printf("\n");

	printf("Msg Encoded after Tempered: [");
	for (size_t i = 0; i < msg->used; i++) {
		printf("%u,",msg->array[i]);
	}
	printf("]\n");
	
	synd = rs_calc_syndromes(msg, ENCODE_SYMBOLS, gf_table);
	printf("synd : ");
	for (size_t i = 0; i < synd->used; i++) {
		printf("%u, ",synd->array[i]);
	}
	printf("\n");
	err_loc = rs_find_error_locator(synd, ENCODE_SYMBOLS, 0, gf_table);
	printf("err_loc : ");
	for (size_t i = 0; i < err_loc->used; i++) {
		printf("%u, ",err_loc->array[i]);
	}
	printf("\n");
	pos = rs_find_errors(reverse_arr(err_loc), msg->used, gf_table);
	printf("err_pos : ");
	for (size_t i = 0; i < pos->used; i++) {
		printf("%u, ",pos->array[i]);
	}
	printf("\n");
	rev_pos = reverse_arr(pos);
	printf("Error positions: [");
	for (size_t i = 0; i < rev_pos->used; i++) {
		printf("%u,", rev_pos->array[i]);
	}
	printf("]\n");
	
	struct Array *err_pos = malloc(sizeof(struct Array));
	initArray(err_pos, ERRORS);
	err_pos->array[0] = 0;
	
	msg = rs_correct_msg(msg, ENCODE_SYMBOLS, err_pos, gf_table);
	
	printf("Msg Corrected: [");
	for (size_t i = 0; i < msg->used; i++) {
		printf("%u,",msg->array[i]);
	}
	printf("]\n");
	
	printf("Msg Corrected: ");
	for (size_t i = 0; i < strlen(my_msg); i++) {
		printf("%c",msg->array[i]);
	}
	printf("\n");

	bool success = true;
	for (size_t i = 0; i < msg_in->used; i++){
		if (msg->array[i] != msg_in->array[i]){
			success = false;
			break;
		}
	}
	if (success){
		printf("Correction success!");
	}
	else{
		printf("Correction failed!");
	}
	printf("\n");

	freeArray(gf_table->gf_exp);
	freeArray(gf_table->gf_log);
	freeArray(msg_in);
	freeArray(msg);
	freeArray(synd);
	freeArray(pos);
	freeArray(rev_pos);

	return 0;
}