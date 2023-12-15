#include "dictionaryEntry.H"

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class stateModel;
Istream& operator>>(Istream&, stateModel&);
Ostream& operator<<(Ostream&, const stateModel&);


class stateModel
{

    private:
        //-
        word name_;

        //-
        word type_;
    
        //-
        word patchName_;


    public:
    
        /// Constructors
    
        //-
        stateModel()
        :
            name_(word::null),
            type_(word::null),
            patchName_(word::null)
        {
        };

        stateModel(const stateModel& sm)
        :
            name_(sm.name_),
            type_(sm.type_),
            patchName_(sm.patchName_)
        {
        };

        //-
        stateModel(Istream& is)
        :
            name_(word::null),
            type_(word::null),
            patchName_(word::null)
        {
            const dictionaryEntry entry(dictionary::null, is);

            name_ = entry.keyword();
            entry.lookup("type") >> type_;
            entry.lookup("patchName") >> patchName_;
        };

        /// Access functions

        //-
        const word& name() const
        {
            return name_;
        }

        //-
        const word& type() const
        {
            return type_;
        };

        //-
        const word& patchName() const
        {
            return patchName_;
        }


        /// Operators
       friend bool operator==
       (
           const stateModel& a,
           const stateModel& b
       )
       {
           return
               (a.name_ == b.name_)
            && (a.type_ == b.type_)
            && (a.patchName_ == b.patchName_);
       }

       friend bool operator!=
       (
            const stateModel& a,
            const stateModel& b
       )
       {
               return !(a == b);
       }


        /// Friend operators

        friend Istream& operator>>(Istream&, stateModel&);
        friend Ostream& operator<<(Ostream&, const stateModel&);
};

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //
/*
Istream& operator>>(Istream& is, specieModel& sm)
{
    const dictionaryEntry entry(dictionary::null, is);

    sm.name_ = entry.keyword();
    entry.lookup("type") >> sm.type_;
    entry.lookup("patchName") >> sm.patchName_;
    
    return is;
}

Ostream& operator<<(Ostream& os, const specieModel& sm)
{
    os  << sm.name_ << nl << token::BEGIN_BLOCK << nl;
    os  << "type " << sm.type_ << token::END_STATEMENT << nl;
    os  << sm.patchName_ << token::END_STATEMENT << nl;
    os  << token::END_BLOCK << nl;
    
    return os;
}
*/
}
