/*
 * aa_aafcommon.cpp -- Common functions used to manipulate AAF
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
 * Copyright (c) 2004 LIRIS (University Claude Bernard Lyon 1)
 * Copyright (c) 2005 Nathan Hurst
 *
 * This file is part of libaffa.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with libaffa; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include "aa.h"
#include <cstdio>
#include <iostream>
#include <cassert>

unsigned AAF::last = 0; // at beginnnig


// Create an AAF from an array of doubles
// ! For debug purposes

AAF:: AAF(double v0, const double * t1, const unsigned * t2, unsigned T) 
    : special(AAF_TYPE_AFFINE)
{
    assert(T > 0);

    length=T;
    cvalue=v0;
    coefficients = new double [length];
    indexes = new unsigned [length];

    for (unsigned i = 0; i < length; i++)
    {
        coefficients[i]=t1[i];
        indexes[i]=t2[i];
    }

    if (indexes[length-1] > last) set_default(indexes[length-1]);

}


// Copy constructor

AAF:: AAF(const AAF &P) 
    : special(P.special), coefficients(NULL), indexes(NULL)
{
    unsigned plength = P.get_length();
    cvalue = P.cvalue;
    length = plength;

    if (plength == 0) return;

    coefficients = new double [plength];
    indexes = new unsigned [plength];

    for (unsigned i = 0; i<plength; i++)
    {
        coefficients[i] = P.coefficients[i];
        indexes[i] = P.indexes[i];
    }

}


// Create an AAF from an interval

AAF:: AAF(interval iv) {
    unsigned en = inclast();
    coefficients = new double [1];
    indexes = new unsigned [1];
  
    if(iv.width() == INFINITY) {
        cvalue = 0;
        length = 1;
        coefficients[0] = INFINITY;
        indexes[0] = en;
        special = AAF_TYPE_INFINITE;
    } else {
        cvalue=(iv.right()+iv.left())/2;
        length = 1;
    
        coefficients[0]=(iv.right()-iv.left())/2;
        indexes[0]=en;
        special = AAF_TYPE_AFFINE;
    }
}


// AAF destructor

AAF::~AAF()
{
    if (length)
    {
        delete [] coefficients;
        delete [] indexes;
    }
}


//  Affectation operator

AAF & AAF::operator = (const AAF & P)
{
    special = P.special;
    unsigned plength = P.get_length();
    
    if (&P!=this)
    {
        if (length != plength)
        {

            if (length)
            {
                delete [] coefficients;
                delete [] indexes;
            }

            if (plength > 0) {
                coefficients = new double [plength];
                indexes = new unsigned [plength];
            }
            else {
                coefficients = NULL;
                indexes = NULL;
            }
        }

        cvalue = P.cvalue;
        length=plength;
        for (unsigned i = 0; i<plength; i++)
        {
            coefficients[i]=P.coefficients[i];
            indexes[i]=P.indexes[i];
        }
    }

    return *this;
}

bool AAF::operator == (const AAF & P) const
{
    if(special!=P.special) return false;
    if(cvalue!=P.cvalue) return false;
    if(length!=P.get_length()) return false;
    for(unsigned i=0; i<length; i++)
    {
        if(coefficients[i]!=P.coefficients[i] ||
           indexes[i]!=P.indexes[i])
            return false;
    }
    return true;
}

bool AAF::operator != (const AAF & P) const
{
    return !(*this==P);
}


namespace aaf_possible_compare
{

bool operator<(const AAF& lhs, double rhs)
{
    return lhs.convert().left() < rhs;
}

bool operator>(const AAF& lhs, double rhs)
{
    return lhs.convert().right() < rhs;
}

bool operator<=(const AAF& lhs, double rhs)
{
    return lhs.convert().left() <= rhs;
}

bool operator>=(const AAF& lhs, double rhs)
{
    return lhs.convert().right() >= rhs;
}

bool operator<(const AAF& lhs, const AAF& rhs)
{
    return lhs - rhs < 0.;
}

bool operator>(const AAF& lhs, const AAF& rhs)
{
    return lhs - rhs > 0.;
}

bool operator<=(const AAF& lhs, const AAF& rhs)
{
    return lhs - rhs <= 0.;
}

bool operator>=(const AAF& lhs, const AAF& rhs)
{
    return lhs - rhs >= 0.;
}


} // end namespace



namespace aaf_absolute_compare
{

bool operator<(const AAF& lhs, double rhs)
{
    return lhs.convert().right() < rhs;
}

bool operator>(const AAF& lhs, double rhs)
{
    return lhs.convert().left() < rhs;
}

bool operator<=(const AAF& lhs, double rhs)
{
    return lhs.convert().right() <= rhs;
}

bool operator>=(const AAF& lhs, double rhs)
{
    return lhs.convert().left() <= rhs;
}

bool operator<(const AAF& lhs, const AAF& rhs)
{
    return lhs - rhs < 0.;
}

bool operator>(const AAF& lhs, const AAF& rhs)
{
    return lhs - rhs > 0.;
}

bool operator<=(const AAF& lhs, const AAF& rhs)
{
    return lhs - rhs <= 0.;
}

bool operator>=(const AAF& lhs, const AAF& rhs)
{
    return lhs - rhs >= 0.;
}

} // end namespace


//xxx bool AAF::operator < (const AAF & P) const
//xxx {
//xxx     return (*this).convert().right() < P.convert().left();
//xxx }
//xxx 
//xxx bool AAF::operator > (const AAF & P) const
//xxx {
//xxx     return (*this).convert().left() > P.convert().right();
//xxx }


// Ostream output of an AAF

std::ostream & operator << (std::ostream & s, const AAF &P)
{

    // s.setf(0, ios_base::floatfield);
    s << "-------------\n";
    s << "Length = " << P.length << "\n";
    s << "v0 = " << P.cvalue << "\n";

    for (unsigned i=0; i < P.length ; i++)
        s << "e" << P.indexes[i] << " = " << P.coefficients[i] << "\n";

    return s;

}


// Print length and coefficients of an AAF to stdout

void AAF::aafprint() const
{
    std::cout << "-------------\n";

    printf("Size = %d\n", length);
    printf("v0 = %f\n", cvalue);

    for (unsigned i=0; i < length ; i++)
        printf("e%d = %f\n", indexes[i], coefficients[i]);
}


// Get the central value x0 of an AAF
// not inline as it isn't used by the lib
// but is useful for applications


double AAF::get_center() const
{
    return cvalue;
}


// Convert an AAF to an interval representation

interval AAF::convert() const
{

    // lower bound == central value of the AAF - the total deviation
    // upper bound == central value of the AAF + the total deviation
    if(is_indeterminate())
        return interval(-INFINITY, INFINITY);
    return interval(cvalue-rad(), cvalue+rad());
}


// Get the total deviation of an AAF
// i.e. the sum of all noise symbols (their abs value)

double AAF::rad() const
{
    double sum=0;

    for (unsigned i=0; i< length; i++) {
        if (coefficients[i] >= 0.0)
            sum+=coefficients[i];
        else
            sum+=-coefficients[i];
    }


    return sum;

}


//xxxx // Added by luk:
//xxxx // Get the upper bound of this AAF. At the same time, record the polarity 
//xxxx // (i.e. +1 or -1) of each noise symbol
//xxxx double AAF::max(AAF& P) const
//xxxx {
//xxxx     P = *this; // copy this AAF
//xxxx     P.cvalue = AAF::last;
//xxxx 
//xxxx     double sum=cvalue;
//xxxx 
//xxxx     for (unsigned i=0; i< length; i++) {
//xxxx         if (coefficients[i] >= 0.0) {
//xxxx             sum+=coefficients[i];
//xxxx             P.coefficients[i] = 1;
//xxxx         }
//xxxx         else {
//xxxx             sum+=-coefficients[i];
//xxxx             P.coefficients[i] = -1;
//xxxx         }
//xxxx     }
//xxxx     return sum;
//xxxx }
//xxxx 
//xxxx double AAF::maximum(std::map<unsigned int,int>& S) const
//xxxx {
//xxxx     double sum = get_center();
//xxxx     
//xxxx     for (unsigned i=0; i< get_length(); i++) {
//xxxx         //assert(coefficients[i] != 0);
//xxxx         if (get_coeff(i) > 0.0) {
//xxxx             sum += get_coeff(i);
//xxxx             S[get_index(i)] = 1;
//xxxx         }
//xxxx         else if (get_coeff(i) < 0.0) {
//xxxx             sum -= get_coeff(i);
//xxxx             S[get_index(i)] = -1;
//xxxx         }
//xxxx         else {
//xxxx             S[get_index(i)] = 0;
//xxxx         }
//xxxx     }
//xxxx     return sum;
//xxxx }
//xxxx 
//xxxx // Added by luk
//xxxx // Evaluate this AAF by using the polarity stored in P
//xxxx // Precondition: this AAF has the same noise symbols as P 
//xxxx double AAF::eval(const AAF& P) const
//xxxx {
//xxxx     //xxx assert(P.cvalue == (double) AAF::last);
//xxxx 
//xxxx     std::map<unsigned, unsigned> S; 
//xxxx     // the indexing may be different so that we use a map for the mapping
//xxxx 
//xxxx     double sum=cvalue;
//xxxx     //xxx assert(P.length >= length);
//xxxx 
//xxxx     // Create the mapping
//xxxx     for (unsigned i=0; i< P.length; i++) {
//xxxx         S[P.indexes[i]] = i;
//xxxx     }
//xxxx     for (unsigned i=0; i< length; i++) {
//xxxx         unsigned idx1 = indexes[i];
//xxxx         unsigned idx2 = 0;
//xxxx         // Note: if corresponding index cannot be found, ignore
//xxxx         if (S.find(idx1) != S.end()) {
//xxxx            idx2 = S[idx1];
//xxxx            //assert(idx1 == idx2);
//xxxx            if (P.coefficients[idx2] >= 0.0) {
//xxxx                sum+=coefficients[i];
//xxxx            }
//xxxx            else {
//xxxx                sum+=-coefficients[i];
//xxxx            }
//xxxx         }
//xxxx     }
//xxxx     return sum;
//xxxx }
//xxxx 
//xxxx double AAF::evaluate(const std::map<unsigned int,int>& S) const
//xxxx {
//xxxx     double sum=get_center();
//xxxx  
//xxxx     for (unsigned i=0; i< get_length(); i++) {
//xxxx         unsigned int idx1 = get_index(i);
//xxxx         // Note: if corresponding index cannot be found, throw an exception
//xxxx         if (S.find(idx1) == S.end()) throw;
//xxxx         int polarity = S.find(idx1)->second;
//xxxx         sum += polarity * get_coeff(i);
//xxxx     }
//xxxx     return sum;
//xxxx }



AAF half_plane(const AAF & P) {
    const double a = P.convert().left(); // [a,b] is our interval
    const double b = P.convert().right();
    
    AAF_TYPE type;
    if(P.special == AAF_TYPE_NAN)
        return P;
    
    if(a > 0)
        return P;
    else if(b < 0) {
        type = AAF_TYPE_NAN;
        return AAF(type);
    }
    else if(b == 0) {
        AAF result(0.);
        result.special = (AAF_TYPE)(AAF_TYPE_AFFINE | AAF_TYPE_NAN);
        return result;
    }
    else if(a <= 0) {
        AAF result;
        result.special = (AAF_TYPE)(AAF_TYPE_AFFINE | AAF_TYPE_NAN);
        unsigned plength = P.get_length();
        if (plength > 0) {
            result.coefficients = new double [plength];
            result.indexes = new unsigned [plength];
        }
        else {
            result.coefficients = NULL;
            result.indexes = NULL;
        }            
        result.cvalue = b/2;;
        result.length = plength;
        
        //xxx double rescale = (b - a)/b;
        for (unsigned i = 0; i<plength; i++)
        {
            result.coefficients[i] = (P.coefficients[i]*b);
            result.indexes[i] = P.indexes[i];
        }
        return result;
    }

    return P;
}


/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/


// vim: filetype=c++:expandtab:shiftwidth=4:tabstop=8:softtabstop=4 :
