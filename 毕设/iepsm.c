/*
 * SMART: string matching algorithms research tool.
 * Copyright (C) 2012  Simone Faro and Thierry Lecroq
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * contact the authors at: faro@dmi.unict.it, thierry.lecroq@univ-rouen.fr
 * download the tool at: http://www.dmi.unict.it/~faro/smart/
 *
 * This is an implementation of the EPSM algorithm
 * in S. Faro and O. M. Kulekci.
 */

#include "include/define.h"
#include "include/main.h"

#include <memory.h>
#include <smmintrin.h>
#include <inttypes.h>

#define HASHSIZE 11

typedef union{
   __m128i* symbol16;
   unsigned char* symbol;
} TEXT;

typedef union{
              __m128i  v;
        unsigned  int  ui[4];
   unsigned short int  us[8];
        unsigned char  uc[16];
}VectorUnion;

typedef struct list
{
    struct list *next;
    int pos;
} LIST;

int search4(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
    __m128i* text = (__m128i*)x;
    __m128i* end  = (__m128i*)(x+16*(textlen/16));
    if ((textlen%16)<7) end--;

    int i,count=0;
    VectorUnion P,Z;
    __m128i a,b,p,z;

    Z.ui[0] = Z.ui[1] = Z.ui[2] = Z.ui[3] = 0;
    z= Z.v;
    P.uc[0] = pattern[0];
    P.uc[1] = pattern[1];
    P.uc[2] = pattern[2];
    P.uc[3] = pattern[3];
    p = P.v;

    text++;// leave the naive check of the first block to the end

    while(text<end)
    {
        //check if P[(m-4) ... (m-1)] matches with T[i*16 ... i*16+3], T[i*16+1 ... i*16+4], .... , T[i*16+7 ... i*16+10]
        a      = _mm_mpsadbw_epu8(*text, p, 0x00);
        b      = _mm_cmpeq_epi16(a,z);
        i      = _mm_movemask_epi8(b);
        count += _mm_popcnt_u32(i);

        a      = _mm_blend_epi16(*text,*(text+1),0x0f);
        b      = _mm_shuffle_epi32(a, _MM_SHUFFLE(1,0,3,2));

        //check if P[(m-4) ... (m-1)] matches with T[i*16+8 ... i*16+11], T[i*16+9 ... i*16+12], .... , T[i*16+15 ... i*16+18]
        a      = _mm_mpsadbw_epu8(b, p, 0x00);
        b      = _mm_cmpeq_epi16(a,z);
        i      = _mm_movemask_epi8(b);
        count += _mm_popcnt_u32(i);
        text++;
    }
    count = count / 2;

    // the ending position of the pattern from the first appropriate position T[patlen-1] to the third position of the next 16-byte block is performed naive
    for(i=3; (i<19) && (i<textlen); i++) // j presents possible end points of the pattern
        if (0==memcmp(pattern,&x[i-3],patlen)) count++;

    // note that at the last iteration of the while loop, we have checked if P ends at positions 0,1,and 2 of the last 16-byte block
    // however, what if the last position of the text is beyond 2?
    // for the possibilities that T ends at positions 3,4,5,6,7,8,9,10,11,12,13,14,and 15, we perform naive checks

    for(i = ((unsigned char*) text)+3-x ; i < textlen ; i++)
    {
        if ( 0 == memcmp(pattern,&x[i-3],4) ) count++;
    }

    return count;
}


int epsm(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
    if (patlen<2)         return search1(pattern, patlen, x, textlen);
    if (patlen==2)        return search2(pattern, patlen, x, textlen);
    if (patlen==3)        return search3(pattern, patlen, x, textlen);
    if (patlen==4)        return search4(pattern, patlen, x, textlen);
    if (patlen>=16)       return search16(pattern, patlen, x, textlen);

    unsigned char* y0;
    int i,j,k,count=0;
    VectorUnion P,zero;
    __m128i res,a,b,z,p;

    __m128i* text = (__m128i*)x;
    __m128i* end  = (__m128i*)(x+16*(textlen/16));
     end--;

    zero.ui[0]=    zero.ui[1]=    zero.ui[2]=    zero.ui[3]=0;  z = zero.v;
    P.uc[0] = pattern[patlen-5];  P.uc[1] = pattern[patlen-4];  P.uc[2] = pattern[patlen-3];  P.uc[3] = pattern[patlen-2];  p = P.v;

    i = (patlen-1) / 16; // i points the first 16-byte block that P may end in
    i++;
    text += i;
    for(k=0; k<(i*16+8)-patlen+1; k++) 	if (0 == memcmp(pattern,x+k,patlen)) count++;

    //the loop checks if pattern ends at the second half of text[i] or at the first half of text[i+1]
    while(text < end)
    {
        //check if P[(m-5) ... (m-2)] matches with T[i*16+4 ... i*16+7], T[i*16+5 ... i*16+8], .... , T[i*16+11 ... i*16+14]
        //note thet this corresponds P ends at T[i*16+8],T[i*16+9],...,T[i*16+15]
        res  = _mm_mpsadbw_epu8(*text, p, 0x04);
        b    = _mm_cmpeq_epi16(res,z);
        j    = _mm_movemask_epi8(b);
        if (j)
        {
            y0 = (unsigned char*)(text) + 9 - patlen;
            if ( (j&3)==3 && !memcmp(pattern,y0,patlen)) count++;
            if ( (j&12)==12 && !memcmp(pattern,y0+1,patlen)) count++;
            if ( (j&48)==48 && !memcmp(pattern,y0+2,patlen)) count++;
            if ( (j&192)==192 && !memcmp(pattern,y0+3,patlen)) count++;
            if ( (j&768)==768 && !memcmp(pattern,y0+4,patlen)) count++;
            if ( (j&3072)==3072 && !memcmp(pattern,y0+5,patlen)) count++;
            if ( (j&12288)==12288 && !memcmp(pattern,y0+6,patlen)) count++;
            if ( (j&49152)==49152 && !memcmp(pattern,y0+7,patlen)) count++;
        }

        a   = _mm_blend_epi16(*text,*(text+1),0x0f);
        b   = _mm_shuffle_epi32(a,_MM_SHUFFLE(1,0,3,2));

        //check if P ends at T[(i+1)*16+8],T[(i+1)*16+9],...,T[(i+1)*16+15]
        res  = _mm_mpsadbw_epu8(b, p, 0x04);
        b    = _mm_cmpeq_epi16(res,z);
        j    = _mm_movemask_epi8(b);

        if (j)
        {
            y0 = (unsigned char*)(text) + 9 + 8 - patlen;
            if ( (j&3)==3 && !memcmp(pattern,y0,patlen)) count++;
            if ( (j&12)==12 && !memcmp(pattern,y0+1,patlen)) count++;
            if ( (j&48)==48 && !memcmp(pattern,y0+2,patlen)) count++;
            if ( (j&192)==192 && !memcmp(pattern,y0+3,patlen)) count++;
            if ( (j&768)==768 && !memcmp(pattern,y0+4,patlen)) count++;
            if ( (j&3072)==3072 && !memcmp(pattern,y0+5,patlen)) count++;
            if ( (j&12288)==12288 && !memcmp(pattern,y0+6,patlen)) count++;
            if ( (j&49152)==49152 && !memcmp(pattern,y0+7,patlen)) count++;
        }
        text++;
    }

    for(k = ((unsigned char*)text)+8-x ; k < textlen ; k++)
    {
        if ( 0 == memcmp(pattern,&x[k-patlen+1],patlen) ) count++;
    }

    return count;
}


