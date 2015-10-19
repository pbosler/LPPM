#ifndef MODEL_LOOKUP_TABLES_H
#define MODEL_LOOKUP_TABLES_H

/*
 FUNCTION PROTOTYPES
*/ 
void GetColorForValue_BlueWhiteRed(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_ColdToHot(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_JetLat(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_MaizeBlue(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_NCLDetail(double val, double &r, double &g, double &b, double &min,double &max);
void GetColorForValue_BlueYellowRed(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_BlueWhiteOrangeRed(double val, double &r, double &g, double &b, double &min, double &max);
void GetColorForValue_RainbowDesaturated(double val, double &r, double &g, double &b, double &min, double &max);

#endif