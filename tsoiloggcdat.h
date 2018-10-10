#if!defined(TSOILOGGCDAT_H)
#define TSOILOGGCDAT_H

class Tsoiloggcdata
{

public:
	Tsoiloggcdata(void);


	int get(ifstream& infile);
	int getdel(FILE* infile);

	float col;
	float row;
	char varname[9];
	int carea;
	double orgcarbon;
	char contnent[9];

private:

	int oggcend;
	long curpos;
	long lagpos;

};

#endif

