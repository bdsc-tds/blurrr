#ifndef BIN_XENIUMASSIGNMENT_HPP
#define BIN_XENIUMASSIGNMENT_HPP

#include <vector>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

template <typename M, typename S, typename C>
class Assignment{
private:
    std::vector<M> assigned_mole;
    std::vector<M> ambi_assigned_mole;

    std::vector<S> assigned_to;
    std::vector<S> ambi_assigned_to;

    std::map<S, C> assign_to_count;

public:
    // constructors
    Assignment() = default; // default constructor without any arguments
    Assignment(const Assignment &) = default; // default copy constructor
    Assignment(Assignment &&) = default; // default move constructor
    Assignment & operator=(const Assignment &) = default; // default copy-assignment constructor
    Assignment & operator=(Assignment &&) = default; // default move-assignment constructor

    Assignment(const M &, const S &);
    Assignment(M &&, S &&);

    Assignment(const M &, const S &, const std::vector<M> &, const std::vector<S> &);
    Assignment(M &&, S &&, std::vector<M> &&, std::vector<S> &&);

    // getters
    const std::vector<M> &get_assigned_mole() const;
    const std::vector<M> &get_ambi_assigned_mole() const;
    const std::vector<S> &get_assigned_to() const;
    const std::vector<S> &get_ambi_assigned_to() const;
    const std::map<S, C> &get_assign_to_count() const;

    std::size_t get_assigned_num() const;
    std::size_t get_ambi_assigned_num() const;
    std::size_t get_assign_to_count_num() const;

    // others
    void add_assigned(const M &, const S &);
    void add_assigned(const std::vector<M> &, const std::vector<S> &);
    void add_ambi_assigned(const std::vector<M> &, const std::vector<S> &);
    void add_assign_to_count(const S &, const C &);
    void add_assign_to_count(const std::map<S, C> &);

    // operator overloading
    friend Assignment<M, S, C> operator+(const Assignment<M, S, C> &v1, const Assignment<M, S, C> &v2) {
        Assignment<M, S, C> other;

        other.add_assigned(v1.get_assigned_mole(), v1.get_assigned_to());
        other.add_assigned(v2.get_assigned_mole(), v2.get_assigned_to());

        other.add_ambi_assigned(v1.get_ambi_assigned_mole(), v1.get_ambi_assigned_to());
        other.add_ambi_assigned(v2.get_ambi_assigned_mole(), v2.get_ambi_assigned_to());

        other.add_assign_to_count(v1.get_assign_to_count());
        other.add_assign_to_count(v2.get_assign_to_count());

        return std::move(other);
    }

    Assignment<M, S, C> &operator+=(const Assignment<M, S, C> &);
};

template <typename M, typename S, typename C>
Assignment<M, S, C>::Assignment(const M &mole, const S &to) {
    this->assigned_mole.emplace_back(mole);
    this->assigned_to.emplace_back(to);
}

template <typename M, typename S, typename C>
Assignment<M, S, C>::Assignment(M &&mole, S &&to) {
    this->assigned_mole.emplace_back(std::move(mole));
    this->assigned_to.emplace_back(std::move(to));
}

template <typename M, typename S, typename C>
Assignment<M, S, C>::Assignment(
    const M &mole,
    const S &to,
    const std::vector<M> &ambi_mole,
    const std::vector<S> &ambi_to
) {
    this->assigned_mole.emplace_back(mole);
    this->assigned_to.emplace_back(to);
    
    this->ambi_assigned_mole.insert(
        this->ambi_assigned_mole.end(),
        ambi_mole.begin(),
        ambi_mole.end()
    );
    this->ambi_assigned_to.insert(
        this->ambi_assigned_to.end(),
        ambi_to.begin(),
        ambi_to.end()
    );
}

template <typename M, typename S, typename C>
Assignment<M, S, C>::Assignment(
    M &&mole,
    S &&to,
    std::vector<M> &&ambi_mole,
    std::vector<S> &&ambi_to
) {
    this->assigned_mole.emplace_back(std::move(mole));
    this->assigned_to.emplace_back(std::move(to));
    
    this->ambi_assigned_mole.insert(
        this->ambi_assigned_mole.end(),
        ambi_mole.begin(),
        ambi_mole.end()
    );
    this->ambi_assigned_to.insert(
        this->ambi_assigned_to.end(),
        ambi_to.begin(),
        ambi_to.end()
    );
}

template <typename M, typename S, typename C>
const std::vector<M> &
Assignment<M, S, C>::get_assigned_mole() const {
    return this->assigned_mole;
}

template <typename M, typename S, typename C>
const std::vector<M> &
Assignment<M, S, C>::get_ambi_assigned_mole() const {
    return this->ambi_assigned_mole;
}

template <typename M, typename S, typename C>
const std::vector<S> &
Assignment<M, S, C>::get_assigned_to() const {
    return this->assigned_to;
}

template <typename M, typename S, typename C>
const std::vector<S> &
Assignment<M, S, C>::get_ambi_assigned_to() const {
    return this->ambi_assigned_to;
}

template <typename M, typename S, typename C>
const std::map<S, C> &
Assignment<M, S, C>::get_assign_to_count() const {
    return this->assign_to_count;
}

template <typename M, typename S, typename C>
std::size_t
Assignment<M, S, C>::get_assigned_num() const {
    if (this->assigned_mole.size() != this->assigned_to.size()) {
        throw std::runtime_error("Internal error!");
    }

    return this->assigned_mole.size();
}

template <typename M, typename S, typename C>
std::size_t
Assignment<M, S, C>::get_ambi_assigned_num() const {
    if (this->ambi_assigned_mole.size() != this->ambi_assigned_to.size()) {
        throw std::runtime_error("Internal error!");
    }

    return this->ambi_assigned_mole.size();
}

template <typename M, typename S, typename C>
std::size_t
Assignment<M, S, C>::get_assign_to_count_num() const {
    return this->assign_to_count.size();
}

template <typename M, typename S, typename C>
void
Assignment<M, S, C>::add_assigned(const M &mole, const S &to) {
    this->assigned_mole.emplace_back(mole);
    this->assigned_to.emplace_back(to);

    const auto it = this->assign_to_count.find(to);
    if (it == this->assign_to_count.end()) {
        this->assign_to_count.emplace(
            std::make_pair(
                to,
                static_cast<C>(1)
            )
        );
    } else {
        (it->second)++;
    }
}

template <typename M, typename S, typename C>
void
Assignment<M, S, C>::add_assigned(const std::vector<M> &mole, const std::vector<S> &to) {
    if (mole.size() != to.size()) {
        throw std::invalid_argument("Unmatched size!");
    }

    auto it_mole = mole.begin(), it_to = to.begin();

    while (it_mole != mole.end() && it_to != to.end()) {
        this->add_assigned(*it_mole, *it_to);

        it_mole++;
        it_to++;
    }
}

template <typename M, typename S, typename C>
void
Assignment<M, S, C>::add_ambi_assigned(
    const std::vector<M> &ambi_mole,
    const std::vector<S> &ambi_to
) {
    this->ambi_assigned_mole.insert(
        this->ambi_assigned_mole.end(),
        ambi_mole.begin(),
        ambi_mole.end()
    );

    this->ambi_assigned_to.insert(
        this->ambi_assigned_to.end(),
        ambi_to.begin(),
        ambi_to.end()
    );
}

template <typename M, typename S, typename C>
void
Assignment<M, S, C>::add_assign_to_count(const S &key, const C &val) {
    const auto __existing_it = this->assign_to_count.find(key);

    if (__existing_it == this->assign_to_count.end()) {
        this->assign_to_count.emplace(
            std::make_pair(
                key,
                val
            )
        );
    } else {
        __existing_it->second += val;
    }
}

template <typename M, typename S, typename C>
void
Assignment<M, S, C>::add_assign_to_count(const std::map<S, C> &count) {
    for (auto it = count.begin(); it != count.end(); it++) {
        this->add_assign_to_count(it->first, it->second);
    }
}

template <typename M, typename S, typename C>
Assignment<M, S, C> &
Assignment<M, S, C>::operator+=(const Assignment<M, S, C> &v) {
    this->add_assigned(v.get_assigned_mole(), v.get_assigned_to());
    this->add_ambi_assigned(v.get_ambi_assigned_mole(), v.get_ambi_assigned_to());
    this->add_assign_to_count(v.get_assign_to_count());

    return *this;
}

#endif // BIN_XENIUMASSIGNMENT_HPP
