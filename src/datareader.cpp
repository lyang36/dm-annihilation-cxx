#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "tipsydefs.h"
#include "datareader.h"

using namespace std;

//read native data
DataReader::DataReader(string path, bool isMask, string maskfile){
    filePath_ = path;
    memBuffer_ = 10000000;
    memCursor_ = 0;
    readCursor_ = 0;
    partNumbers_ = 0;
    memParts_ = 0;
    testnum_ = -1;
    buffer_ = NULL;
    
    isNative_ = true;
    isMask_ = isMask;
    mask_file_name = maskfile;
    
    isSingleFile_ = false;
    
    
}

//read XDR data
DataReader::DataReader(string basedir, string basename, bool isMask, string maskfile, bool isSingleFile){
    basedir_ = basedir;
    basename_ = basename;
    
    memBuffer_ = 10000000;
    memCursor_ = 0;
    readCursor_ = 0;
    partNumbers_ = 0;
    memParts_ = 0;
    testnum_ = -1;
    buffer_ = NULL;
    
    isNative_ = false;
    isMask_ = isMask;
    isSingleFile_ = isSingleFile;
    mask_file_name = maskfile;
    
    particle_file_name = basedir_ + "/" + basename_;
    density_file_name = particle_file_name + ".den32";
    hsmooth_file_name = particle_file_name + ".hsm32";
    sigmav_file_name = particle_file_name + ".sigma32";
}

void DataReader::setTest(int tn){
    testnum_ = tn;
}

//open the file and load buffer
bool DataReader::open(){
    readCursor_ = 0;
    if(buffer_ != NULL){
        //free buffer
        delete(buffer_);
    }
    buffer_ = NULL;
    memCursor_ = 0;
    memParts_ = 0;
    buffer_ = new DMParticle[memBuffer_];
    
    
    if(isNative_){
        dataInputFile_.open(filePath_.c_str(), ios::binary);
        
        if(!dataInputFile_.good()){
            printf("Datafile error: %s !\n", filePath_.c_str());
            exit(1);
        }
        dataInputFile_.read((char*)&partNumbers_, sizeof(int));
        cout << "Particles: " << partNumbers_ << endl;
    }else{
        //open particle file
        char fi[255];
        strcpy(fi, particle_file_name.c_str());
        fp_=fopen(fi,"r");
        if(!fp_) {
            fprintf(stderr,"Problem opening file: %s\n",particle_file_name.c_str());
            exit(1);
        }
        xdrstdio_create(&xdrs_,fp_, XDR_DECODE);
        //read header
        int status = xdr_header(&xdrs_, &tips_);
        partNumbers_ = tips_.ndark;
        
        cout << "Particles: " << partNumbers_ << endl;
        
        int np;
        if(!isSingleFile_){
            densStream_.open( density_file_name.c_str(), ios::in | ios::binary);
            densStream_.read( reinterpret_cast<char*>( &np ), sizeof np );
            if(np != partNumbers_){
                fprintf(stderr,"Particle numbers in density not match.\n");
                exit(1);
            }
            
            //open hsmooth file
            hsmoothStream_.open( hsmooth_file_name.c_str(), ios::in | ios::binary);
            hsmoothStream_.read( reinterpret_cast<char*>( &np ), sizeof np );
            if(np != partNumbers_){
                fprintf(stderr,"Particle numbers in hsmooth not match.\n");
                exit(1);
            }
            
            sigmavStream_.open( sigmav_file_name.c_str(), ios::in | ios::binary);
            sigmavStream_.read( reinterpret_cast<char*>( &np ), sizeof np);
            if(np != partNumbers_){
                fprintf(stderr,"Particle numbers in sigmav not match.\n");
                exit(1);
            }
            
        }
    }
    
    if(isMask_){
        int np;
        maskStream_.open( mask_file_name.c_str(), ios::in | ios::binary);
        maskStream_.read( reinterpret_cast<char*>( &np ), sizeof np );
        if(np != partNumbers_){
            fprintf(stderr,"Particle numbers in mask file not match.\n");
            exit(1);
        }
    }
    
    cout << "Load to buffer ... " << endl;
    //loadBuffer();
    if(testnum_ != -1){
        if((int)partNumbers_ > testnum_){
            partNumbers_ = testnum_;
        }
    }
    loadBuffer();
    
    
    if(isNative_){
        return dataInputFile_.is_open();
    }else{
        return (fp_ != 0);
    }
}

void DataReader::loadBuffer(){
    int resParts = partNumbers_ - (readCursor_ - memCursor_ + memParts_);
    //cout << "resparts: " << resParts << endl;
    if(resParts > (int)memBuffer_){
        //setup new readCursor_
        readCursor_ = readCursor_ - memCursor_ + memParts_;
        memCursor_ = 0;
        memParts_ = memBuffer_;
    }else{
        //setup new readCursor_
        readCursor_ = readCursor_ - memCursor_ + memParts_;
        memCursor_ = 0;
        memParts_ = resParts;
    }
    
    //load buffer
    if(isNative_){
        dataInputFile_.read((char *) buffer_, sizeof(DMParticle) * memParts_);
    }else{
        int i;
        Pdm dp;
        float dens = 0;
        float hsmooth = 0;
        float sigmav = 0;
        int status;
        for(i = 0; i < (memParts_); i++){
            status = xdr_dark(&xdrs_, &(dp));
            if (status != 1) {
                fprintf(stderr,"Error reading dark particle from input file.\n");
                exit(1);
            }
            if(!isSingleFile_){
                densStream_.read((char *) &dens, sizeof(float));
                hsmoothStream_.read((char *) &hsmooth, sizeof(float));
                sigmavStream_.read((char *) &sigmav, sizeof(float));
            }
            buffer_[i].setPdm(dp);
            buffer_[i].dens = dens;
            buffer_[i].hsmooth = hsmooth;
            buffer_[i].sigmav = sigmav;
        }
    }
    //apply the mask file
    if(isMask_){
        char * mask = new char[memParts_];
        maskStream_.read(mask, sizeof(char) * memParts_);
        for(int i = 0; i < memParts_; i++){
            if(mask[i] == 0){
                buffer_[i].dens = -1;
            }
        }
    }
}

bool DataReader::isOpen(){
    if(isNative_){
        return dataInputFile_.is_open();
    }else{
        return (fp_ != 0);
    }
}

void DataReader::setBuf(unsigned int mem){
    if(buffer_ != NULL){
        //free buffer_
        delete(buffer_);
    }
    buffer_ = NULL;
    memCursor_ = 0;
    memParts_ = 0;
    memBuffer_ = mem;
}

bool DataReader::hasNext(){
    if(readCursor_ < partNumbers_)
        return true;
    else
        return false;
}

bool DataReader::readParticle(DMParticle * part){
    if(memCursor_ == memBuffer_){
        //cout << "load buffer\n";
        loadBuffer();
    }
    //printf("%d\n", readCursor_);
    if(hasNext()){
        //(*part) = buffer_[memCursor_];
        memcpy((char *)part, buffer_ + memCursor_, sizeof(DMParticle));
        readCursor_ ++; //move the cursor
        memCursor_ ++;  //move the cursor
    }else{
        return false;
    }
    return true;
}

void DataReader::move2bufEnd(){
    readCursor_ += memParts_ - memCursor_;
    memCursor_ = memParts_;
}

void DataReader::close(){
    if(buffer_ != NULL)
        delete[](buffer_);
    buffer_ = NULL;
    if(isNative_){
        dataInputFile_.close();
    }else{
        xdr_destroy(&xdrs_);
        fclose(fp_);
        
        if(!isSingleFile_){
            densStream_.close();
            hsmoothStream_.close();
            sigmavStream_.close();
        }
    }
    if(isMask_){
        maskStream_.close();
    }
}


DMParticle * DataReader::getBuf(){
    return buffer_;
}

unsigned int DataReader::getBufSize(){
    return memParts_;
}
