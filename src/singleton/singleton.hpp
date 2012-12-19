#ifndef _SINGLETON_HPP_
#define _SINGLETON_HPP_

template <typename T>
class SingletonHolder {
    public:
        static T& instance() {
            static T obj;
            return obj;
        }
    private:
        SingletonHolder();
        SingletonHolder(const SingletonHolder& );
        SingletonHolder& operator=(const SingletonHolder& );
        ~SingletonHolder();
};


#endif // _SINGLETON_HPP_
