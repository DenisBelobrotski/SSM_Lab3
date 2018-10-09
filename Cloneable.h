#ifndef LAB3_CLONEABLE_H
#define LAB3_CLONEABLE_H

class Cloneable
{
public:
    virtual Cloneable* clone() const = 0;
};

#endif //LAB3_CLONEABLE_H
