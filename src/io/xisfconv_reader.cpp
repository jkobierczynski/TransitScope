// SPDX-License-Identifier: GPL-3.0-or-later
// libxisfconv backend (FITS, XISF, ASDF, DNG through one API).
//
// Mapping onto ImageReader:
//   hdu index  -> xisfconv image index (images only: tables and other skipped HDUs are not counted)
//   header()   -> xisfconv_image_keywords: the non-structural cards (SIMPLE, BITPIX, NAXISn, BZERO,
//                 BSCALE, ... are not part of it; fitscompare reports those separately)
//   pixels()   -> samples as stored (integer types are physical values after BZERO/BSCALE, as the
//                 library maps them), converted to float here; rows in FITS order (bottom-up), so
//                 FITS files match CFITSIO and XISF files come out in the orientation their WCS
//                 keywords describe.
#include <xisfconv.h>

#include <algorithm>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <string>

#include "transit/image_reader.h"

namespace transit {

HeaderCard::Type classify_fits_value(const std::string& value);  // reader_factory.cpp

namespace {

class XisfconvReader final : public ImageReader {
public:
    XisfconvReader() : ctx_(xisfconv_context_new()) {
        if (!ctx_) throw std::runtime_error("xisfconv: cannot create context");
    }
    ~XisfconvReader() override {
        xisfconv_close(file_);
        xisfconv_context_free(ctx_);
    }

    void open(const std::string& path) override {
        xisfconv_close(file_);
        file_ = nullptr;
        check(xisfconv_open(ctx_, path.c_str(), &file_), "open " + path);
    }

    int hdu_count() const override { return static_cast<int>(xisfconv_image_count(need())); }

    std::vector<HeaderCard> header(int hdu) const override {
        const xisfconv_keywords* kw = nullptr;
        check(xisfconv_image_keywords(need(), static_cast<size_t>(hdu), &kw), "keywords");
        std::vector<HeaderCard> cards;
        const size_t n = xisfconv_keywords_count(kw);
        cards.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            const char *name = nullptr, *value = nullptr, *comment = nullptr;
            check(xisfconv_keywords_get(kw, i, &name, &value, &comment), "keyword");
            HeaderCard c{name ? name : "", value ? value : "", comment ? comment : "", HeaderCard::Type::Undefined};
            c.type = classify_fits_value(c.raw_value);
            cards.push_back(std::move(c));
        }
        return cards;
    }

    std::vector<long> dims(int hdu) const override {
        const auto info = image_info(hdu);
        std::vector<long> d{static_cast<long>(info.width), static_cast<long>(info.height)};
        if (info.channels > 1) d.push_back(static_cast<long>(info.channels));
        return d;
    }

    // Read in bands of rows into a small reused buffer and convert straight into the Image: the
    // library then reads the data where it is in the file instead of loading and copying the
    // whole image twice (about 3x faster on uncompressed FITS, same values).
    Image pixels(int hdu) const override {
        const size_t idx = static_cast<size_t>(hdu);
        xisfconv_read_options opt;
        xisfconv_read_options_init(&opt, sizeof opt);
        opt.sample_format = XISFCONV_SAMPLE_AS_STORED;
        opt.row_order = XISFCONV_ROWS_BOTTOM_UP;
        check(xisfconv_read_rows(need(), idx, &opt, 0, 0, nullptr, 0), "read_rows");  // decides the sample format
        const auto info = image_info(hdu);
        if (!info.convertible)
            throw std::runtime_error(std::string("xisfconv: image not readable: ") +
                                     xisfconv_image_unsupported_reason(need(), idx));
        const size_t w = info.width, h = info.height, ch = std::max<size_t>(1, info.channels);
        const size_t sb = xisfconv_sample_size(info.sample_format);
        if (!sb) throw std::runtime_error("xisfconv: unknown sample format");
        constexpr size_t kBand = 64;
        band_.resize(kBand * w * ch * sb);

        Image img(static_cast<int>(w), static_cast<int>(h));
        for (size_t y = 0; y < h; y += kBand) {
            const size_t rows = std::min(kBand, h - y);
            check(xisfconv_read_rows(need(), idx, &opt, y, rows, band_.data(), band_.size()), "read_rows");
            float* dst = img.px.data() + y * w;
            const size_t n = rows * w;  // channel 0 comes first in the planar band
            auto conv = [&](auto* src) { for (size_t i = 0; i < n; ++i) dst[i] = static_cast<float>(src[i]); };
            const unsigned char* buf = band_.data();
            switch (info.sample_format) {
                case XISFCONV_SAMPLE_UINT8:   conv(reinterpret_cast<const uint8_t*>(buf)); break;
                case XISFCONV_SAMPLE_UINT16:  conv(reinterpret_cast<const uint16_t*>(buf)); break;
                case XISFCONV_SAMPLE_UINT32:  conv(reinterpret_cast<const uint32_t*>(buf)); break;
                case XISFCONV_SAMPLE_UINT64:  conv(reinterpret_cast<const uint64_t*>(buf)); break;
                case XISFCONV_SAMPLE_FLOAT32: conv(reinterpret_cast<const float*>(buf)); break;
                case XISFCONV_SAMPLE_FLOAT64: conv(reinterpret_cast<const double*>(buf)); break;
                default: throw std::runtime_error("xisfconv: unknown sample format");
            }
        }
        return img;
    }

    std::string backend() const override { return "xisfconv"; }

private:
    xisfconv_image_info image_info(int hdu) const {
        xisfconv_image_info info;
        xisfconv_image_info_init(&info, sizeof info);
        check(xisfconv_image_info_get(need(), static_cast<size_t>(hdu), &info), "image_info");
        return info;
    }
    xisfconv_file* need() const {
        if (!file_) throw std::runtime_error("xisfconv: no file open");
        return file_;
    }
    void check(xisfconv_status st, const std::string& what) const {
        if (st == XISFCONV_OK) return;
        const char* msg = xisfconv_error_message(ctx_);
        throw std::runtime_error("xisfconv " + what + ": " + xisfconv_status_text(st) +
                                 (msg && *msg ? std::string(" (") + msg + ")" : ""));
    }

    xisfconv_context* ctx_ = nullptr;
    xisfconv_file* file_ = nullptr;
    mutable std::vector<unsigned char> band_;  // per reader, so per thread
};

}  // namespace

std::unique_ptr<ImageReader> make_xisfconv_reader() { return std::make_unique<XisfconvReader>(); }

}  // namespace transit
