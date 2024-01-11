#ifndef BINXENIUM_SUBSPOTS_HPP
#define BINXENIUM_SUBSPOTS_HPP

template<typename T>
class Subspot {
private:
    T idx;
    T id;
    T third_pt_id;

public:
    // constructors
    Subspot() = delete; // disable the default constructor without any arguments

    Subspot(const Subspot &) = default; // default copy constructor
    Subspot(Subspot &&) = default; // default move constructor
    Subspot & operator=(const Subspot &) = default; // default copy-assignment constructor
    Subspot & operator=(Subspot &&) = default; // default move-assignment constructor

    Subspot(T, T, T);

    // getters
    T get_idx() const;
    T get_id() const;
    T get_third_pt_id() const;
};

template<typename T>
Subspot<T>::Subspot(T idx, T id, T third_pt_id): idx(idx), id(id), third_pt_id(third_pt_id) {}

template<typename T>
T
Subspot<T>::get_idx() const {
    return this->idx;
}

template<typename T>
T
Subspot<T>::get_id() const {
    return this->id;
}

template<typename T>
T
Subspot<T>::get_third_pt_id() const {
    return this->third_pt_id;
}

#endif // BINXENIUM_SUBSPOTS_HPP
