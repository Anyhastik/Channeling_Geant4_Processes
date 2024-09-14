// #pragma once
#ifndef PHYS_CONST_h
#define PHYS_CONST_h 


#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#define ld1 long double
//#define d double
#include "globals.hh"
using namespace std;

constexpr ld1 c0 = 299792458;  /// �������� ����� m/s;
#define echarge (4.8032 * pow(10, -10) * sqrt(1 / (1.60217733 * pow(10, -6)) * pow(10, -2))) /// ����� ���������
constexpr int nmax = 20;  /// ����� ����� � ���������� ��������� ;
#define hbar ((6.582) * (pow(10, -22)))  /// ���������� ������ MeV * s 
#define m0 ((0.5109996) * (pow(10, 6)))  /// ev/c^2
constexpr ld1 e0 = 0.5109996;          /// ������� ����� MeV
#define e (800)                          /// ������� ��������� MeV 
#define mom ((sqrt(pow(e, 2)) - (pow(e0, 2)))) /// MeV
#define gama ((e) / (e0)) /// �����-������
#define coeff ((1.123443973421022) * (pow(10, 28))) /// coeff = - pow(c0, 2) * pow(10, 20) / (gama * m0);
constexpr ld1 tetta = 0;  /// ���� ����� ������� � ��������;
constexpr int nbeam = 100;  /// ����� ����� �����;
constexpr int n_steps = 3000;  /// ����� ����� ��� ������� ��������� ��������
const G4String crystaltype = "w";  /// ��� ���������
constexpr ld1 crystaltickness = 40;  /// ������� ���������
/// ������� �������
constexpr int k_orient = 1;
constexpr int l_orient = 1;
constexpr int n_orient = 1;

/**
* ����� �������� ������������ ��� ������� ����� ����������
* � ����������� �� ���� ��������� crystaltype
*/
class GPhysConst
{
public:
	GPhysConst() {}
	virtual ~GPhysConst() {}
	
static tuple<ld1, array<ld1, 2*nmax+1>> cpotss(G4int sigen, G4int k0, G4int l0, G4int n0, G4String crystaltype);
};

#endif  /* PHYS_CONST*/