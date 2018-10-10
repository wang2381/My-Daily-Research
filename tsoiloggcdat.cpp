

/* **********************************************

TSOILOGGCDAT.CPP - object to read the structure of soil organic carbon data for pixels from the previous simulations

************************************************ */

#if!defined(TSOILOGGCDAT_H)
#include"tsoiloggcdat.h"
#endif

/* *************************************************** */

Tsoiloggcdata::Tsoiloggcdata(void)
{
	oggcend = 1;
	lagpos = -99;
	curpos = 0;

};

/* *********************************

public functions

*********************************** */

int Tsoiloggcdata::get(ifstream& infile)
{

	lagpos = infile.tellg();

	infile >> col >> row;
	infile >> varname;
	infile >> carea;
	infile >> orgcarbon;
	infile >> contnent;

infile.seekg(0,ios::cur);
curpos = infile.tellg();

if (curpos < (lagpos+10)) {oggcend = -1;}
return oggcend;

}


int Tsoiloggcdata::getdel(FILE* infile)
{
	oggcend = fscanf(infile, "%f,%f, %s ,%d,%lf, %s",
		     &col,&row,varname,&carea,&orgcarbon,contnent);

	return oggcend;

};

