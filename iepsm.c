#include <x86intrin.h>
#include <memory.h>
#include <smmintrin.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define DNASIZE 2048
typedef struct node
{
    struct node *next;
    int pos;
    unsigned long long val;
} NODE;
/*for genome*/
int iepsm1(unsigned char* p, int m, unsigned char* t, int n ,int alplen){
    NODE* nodes[DNASIZE]; //11 bit hash is gives the best result according to our tests, no shorter no longer
    NODE *t1;
    unsigned int i,filter,shift = m-alplen ;
    unsigned long long crc, seed= 123456789, mask,mask2=0xffffffffffffffff,mask3;
    unsigned long long* ptr64;
    unsigned char* endPtr;
    unsigned char* ptr8;
    unsigned char parr[16];
    int count=0;
    unsigned long long _val;
    const int ALP = 8;
    mask = DNASIZE-1;
    if (m<8){
    	return -1;
    }
    if (m<16){
    	for(i=0;i<m;i++){
    		parr[i] = p[i];
    	}
    	for(i=m;i<16;i++){
    		parr[i] = 0;
    	}
    	p = parr;
    }
    if(alplen==5){
    	mask2 = 0x000000ffffffffff;
    	mask3 = 0xffffff0000000000;
    }else if(alplen == 6){
    	mask2 = 0x0000ffffffffffff;
    	mask3 = 0xffff000000000000;
    }else if(alplen == 7){
    	mask2 = 0x00ffffffffffffff;
    	mask3 = 0xff00000000000000;
    }else{
    	mask2 = 0xffffffffffffffff;
    }
    memset(nodes,0,sizeof(NODE*)*DNASIZE);
    for(i=0;i<m-alplen;i++){
        filter = (unsigned int)(_mm_crc32_u64(seed,(*(unsigned long long*)(&p[i])&mask2)) & mask);
        if(nodes[filter]==0){
        	nodes[filter] = (NODE*)malloc(sizeof(NODE));
            if (nodes[filter]){
	            nodes[filter]->next = 0;
        	    nodes[filter]->pos  = i;
        	    nodes[filter]->val = _mm_crc32_u64(seed, (*(unsigned long long*)(&p[i])&mask2));
            }
        }else{
        	t1 = nodes[filter];
			while (t1->next != 0){
				t1 = t1->next;
			}
			t1->next = (NODE*) malloc(sizeof(NODE));
			if (t1->next) {
				t1 = t1->next;
				t1->next = 0;
				t1->pos = i;
				t1->val = _mm_crc32_u64(seed, (*(unsigned long long*)(&p[i])&mask2));
			}
        }
    }
    endPtr = &t[n-alplen];
    ptr8 = &t[m-alplen];
    while(ptr8 < endPtr){
    	_val = _mm_crc32_u64(seed, (*(unsigned long long*)(ptr8)&mask2));
    	filter = (unsigned int) ( _val & mask);
        if (nodes[filter]){
            t1 = nodes[filter];
            while(t1){
            	if (_val==t1->val){
            		if (memcmp(p,ptr8 - t1->pos,m) == 0){
            		   	count++;
            		}
            	}

                t1=t1->next;
            }
        }
        ptr8+=shift;
    }
    return count;
}
/*for acid*/
int iepsm2(unsigned char* p, int m, unsigned char* t, int n ,int alplen){
    NODE* nodes[DNASIZE]; //11 bit hash is gives the best result according to our tests, no shorter no longer
    NODE *t1;
    alplen=4;
    unsigned int i,filter,shift = m-alplen ;
    unsigned long long crc, seed= 123456789, mask,mask2=0xffffffffffffffff,mask3;
    unsigned long long* ptr64;
    unsigned int* ptr32;
    unsigned char* endPtr;
    unsigned char* ptr8;
    int count=0;
    unsigned long long _val;
    const int ALP = 8;
    mask = DNASIZE-1;
    if (m<8){
    	return -1;
    }
    memset(nodes,0,sizeof(NODE*)*DNASIZE);

    for(i=0;i<m-alplen;i++){
    	filter = (_mm_crc32_u32(seed,(*(unsigned int*)&p[i])) & mask);
        if(nodes[filter]==0){
        	nodes[filter] = (NODE*)malloc(sizeof(NODE));
            if (nodes[filter]){
	            nodes[filter]->next = 0;
        	    nodes[filter]->pos  = i;
        	    nodes[filter]->val = _mm_crc32_u32(seed, *(unsigned int*)(&p[i]));
            }
        }else{
        	t1 = nodes[filter];
			while (t1->next != 0){
				t1 = t1->next;
			}
			t1->next = (NODE*) malloc(sizeof(NODE));
			if (t1->next) {
				t1 = t1->next;
				t1->next = 0;
				t1->pos = i;
				t1->val = _mm_crc32_u32(seed, *(unsigned int*)&p[i]);
			}
        }
    }
    endPtr = &t[n-alplen];
    ptr8 = &t[m-alplen];

    while(ptr8 < endPtr){
    	_val = _mm_crc32_u32(seed, (*(unsigned int*)(ptr8)));
    	filter = ( _val & mask);
        if (nodes[filter]){
            t1 = nodes[filter];
            while(t1){
            	if (_val==t1->val){
            		if (memcmp(p,ptr8 - t1->pos,m) == 0){
            		   	count++;
            		}
            	}

                t1=t1->next;
            }
        }
        ptr8+=shift;
    }
    return count;
}
