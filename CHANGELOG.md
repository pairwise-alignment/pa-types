# CHANGELOG

<!-- next-header -->

## git

## 1.5.0
- Bugfix: Do not panic on case-only substitutions in `to_char_pairs`.

## 1.4.0
- Add `CigarOp::edit_cost` that returns 0 for `Match` and 1 for `Sub`/`Del`/`Ins`.
- Add `Cigar::pop_op` that pops the last operation.
- Improve `CigarOpChars` docs.

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
