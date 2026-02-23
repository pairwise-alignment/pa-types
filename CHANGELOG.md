# CHANGELOG

<!-- next-header -->

## git

## 1.3.0
- Improve docs around the difference between cigar insert (extra in pattern/query) and delete
  (extra in text/reference).
- Improve variable names from `a` and `b` to `text` and `pattern`.
- `Pos` is consistently `(text pos, pattern pos)`.
- Flip `Cigar::to_char_pairs` argument names to be `(text, pattern)`.
  Functionally changes nothing, but now the SAM spec for Cigar insert/delete is followed.

## 1.2.0
- Add `cigar.to_char_pairs(pattern, text)` to iterate corresponding characters.
- Add `cigar.clear()` for `cigar.ops.clear()`.

## 1.1.0
- #2: `cigar.to_string()` explicitly writes counts of 1.

## 1.0.0
- Initial release
