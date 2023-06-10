// Copyright 2007-2010 Baptiste Lepilleur and The JsonCpp Authors
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

// included by json_value.cpp

namespace Json
{

  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // class ValueIteratorBase
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////

  inline ValueIteratorBase::ValueIteratorBase() : current_() {}

  inline ValueIteratorBase::ValueIteratorBase(
      const Value::ObjectValues::iterator &current)
      : current_(current), isNull_(false) {}

  inline Value &ValueIteratorBase::deref() { return current_->second; }
  inline const Value &ValueIteratorBase::deref() const { return current_->second; }

  inline void ValueIteratorBase::increment() { ++current_; }

  inline void ValueIteratorBase::decrement() { --current_; }

  inline ValueIteratorBase::difference_type
  ValueIteratorBase::computeDistance(const SelfType &other) const
  {
    // Iterator for null value are initialized using the default
    // constructor, which initialize current_ to the default
    // std::map::iterator. As begin() and end() are two instance
    // of the default std::map::iterator, they can not be compared.
    // To allow this, we handle this comparison specifically.
    if (isNull_ && other.isNull_)
    {
      return 0;
    }

    // Usage of std::distance is not portable (does not compile with Sun Studio 12
    // RogueWave STL,
    // which is the one used by default).
    // Using a portable hand-made version for non random iterator instead:
    //   return difference_type( std::distance( current_, other.current_ ) );
    difference_type myDistance = 0;
    for (Value::ObjectValues::iterator it = current_; it != other.current_;
         ++it)
    {
      ++myDistance;
    }
    return myDistance;
  }

  inline bool ValueIteratorBase::isEqual(const SelfType &other) const
  {
    if (isNull_)
    {
      return other.isNull_;
    }
    return current_ == other.current_;
  }

  inline void ValueIteratorBase::copy(const SelfType &other)
  {
    current_ = other.current_;
    isNull_ = other.isNull_;
  }

  inline Value ValueIteratorBase::key() const
  {
    const Value::CZString czstring = (*current_).first;
    if (czstring.data())
    {
      if (czstring.isStaticString())
        return Value(StaticString(czstring.data()));
      return Value(czstring.data(), czstring.data() + czstring.length());
    }
    return Value(czstring.index());
  }

  inline UInt ValueIteratorBase::index() const
  {
    const Value::CZString czstring = (*current_).first;
    if (!czstring.data())
      return czstring.index();
    return Value::UInt(-1);
  }

  inline String ValueIteratorBase::name() const
  {
    char const *keey;
    char const *end;
    keey = memberName(&end);
    if (!keey)
      return String();
    return String(keey, end);
  }

  inline char const *ValueIteratorBase::memberName() const
  {
    const char *cname = (*current_).first.data();
    return cname ? cname : "";
  }

  inline char const *ValueIteratorBase::memberName(char const **end) const
  {
    const char *cname = (*current_).first.data();
    if (!cname)
    {
      *end = nullptr;
      return nullptr;
    }
    *end = cname + (*current_).first.length();
    return cname;
  }

  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // class ValueConstIterator
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////

  inline ValueConstIterator::ValueConstIterator() = default;

  inline ValueConstIterator::ValueConstIterator(
      const Value::ObjectValues::iterator &current)
      : ValueIteratorBase(current) {}

  inline ValueConstIterator::ValueConstIterator(ValueIterator const &other)
      : ValueIteratorBase(other) {}

  inline ValueConstIterator &ValueConstIterator::
  operator=(const ValueIteratorBase &other)
  {
    copy(other);
    return *this;
  }

  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // class ValueIterator
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////

  inline ValueIterator::ValueIterator() = default;

  inline ValueIterator::ValueIterator(const Value::ObjectValues::iterator &current)
      : ValueIteratorBase(current) {}

  inline ValueIterator::ValueIterator(const ValueConstIterator &other)
      : ValueIteratorBase(other)
  {
    throwRuntimeError("ConstIterator to Iterator should never be allowed.");
  }

  inline ValueIterator::ValueIterator(const ValueIterator &other) = default;

  inline ValueIterator &ValueIterator::operator=(const SelfType &other)
  {
    copy(other);
    return *this;
  }

} // namespace Json
