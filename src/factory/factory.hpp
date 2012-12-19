#ifndef _FACTORY_HPP_
#define _FACTORY_HPP_

#include <map>

template <typename T>
class Factory {
    public:
        Factory() {}

        typedef T* (*Creator)();

        bool add(const std::string& name, Creator creator) {
            return creators.insert(std::make_pair(name, creator)).second;
        }

        T* create(const std::string& name) {
            typename Creators::const_iterator p = creators.find(name);
            if (p == creators.end()) {
                // TODO: raise exception
            }
            return (p->second)();
        }

    private:
        typedef std::map<std::string, Creator> Creators;
        Creators creators;
};

#endif // _FACTORY_HPP_
